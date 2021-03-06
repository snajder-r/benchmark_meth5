import numpy as np
from benchmark_pycometh.config import module_config
from meth5.meth5 import MetH5File, ChromosomeContainer
from meth5.sparse_matrix import SparseMethylationMatrixContainer


def sample_fun(sample, met_matrix, m5):
    haplotype_map = {int(k): v for k, v in m5.h5_fp["reads/read_groups"]["haplotype"].attrs.items()}
    ret = [haplotype_map.get(x, "none") for x in met_matrix.read_samples]
    return np.array([f"{sample} ({hp})" if hp != "none" else sample for hp in ret])


class MetMatrixLoader:
    def __init__(self):
        self.h5 = {
            sample: MetH5File(module_config.meth5_template_file.format(sample=sample), "r")
            for sample in module_config.samples
        }
    
    def __del__(self):
        for h5file in self.h5.values():
            try:
                h5file.close()
            except:
                pass
    
    def get_merged_matrix(
        self, chrom, start, end, samples=module_config.samples, must_overlap_position=None, sample_name_fun=sample_fun
    ):
        merged_matrix = None
        for sample in module_config.samples:
            if sample not in samples:
                continue
            chrom_container: ChromosomeContainer = self.h5[sample][chrom]
            met_matrix = chrom_container.get_values_in_range(start, end).to_sparse_methylation_matrix(
                read_groups_key="haplotype", read_read_names=False
            )
            
            if must_overlap_position is not None:
                has_over = (
                    np.array(
                        (met_matrix.met_matrix[:, met_matrix.genomic_coord > must_overlap_position] != 0).sum(axis=1)
                    ).flatten()
                    > 0
                )
                has_under = (
                    np.array(
                        (met_matrix.met_matrix[:, met_matrix.genomic_coord < must_overlap_position] != 0).sum(axis=1)
                    ).flatten()
                    > 0
                )
                read_idx = has_over & has_under
                met_matrix = met_matrix.get_submatrix_from_read_mask(read_idx)
                if met_matrix.met_matrix.shape[1] == 0:
                    continue
            
            met_matrix.read_samples = sample_name_fun(sample, met_matrix, self.h5[sample])
            if merged_matrix is None:
                merged_matrix = met_matrix
            else:
                try:
                    merged_matrix = merged_matrix.merge(met_matrix, sample_names="keep")
                except:
                    print("Nothing to merge in")
                    pass
        return merged_matrix
