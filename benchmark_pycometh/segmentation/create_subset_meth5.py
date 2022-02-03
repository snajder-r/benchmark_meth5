import random
from typing import Dict
import tqdm
from meth5.meth5 import MetH5File
import pandas as pd
import numpy as np

from benchmark_pycometh.config import module_config
from benchmark_pycometh.bsseq.create_mock_bs_from_nanopore import MetH5ToBedGraph

ori_m5 = MetH5File(module_config.meth5_template_file.format(sample="HG003"), "r")

values_container = ori_m5["21"].get_all_values()
llrs = values_container.get_llrs()
ranges = values_container.get_ranges()
read_ids = values_container.get_read_ids()

with MetH5File(module_config.subset_m5_file.format(sample="HG003"), "w") as out_m5:
    unique_read_ids = set(read_ids)
    for i in range(100):
        subset_unique_reads = set(random.sample(unique_read_ids, len(unique_read_ids) // 2))
        idx = [r in subset_unique_reads for r in tqdm.tqdm(read_ids)]
        subset_llrs = llrs[idx]
        subset_ranges = ranges[idx]
        subset_read_ids = read_ids[idx]
        
        out_m5.add_to_h5_file(
            pd.DataFrame(
                dict(
                    chromosome=f"21_{i}",
                    log_lik_ratio=subset_llrs,
                    start=subset_ranges[:, 0],
                    end=subset_ranges[:, 1],
                    read_id=subset_read_ids,
                )
            ), postpone_sorting_until_close=True
        )
    ori_m5.h5_fp.copy(ori_m5.h5_fp["reads"], out_m5.h5_fp)


    
converter = MetH5ToBedGraph(module_config.subset_m5_file.format(sample="HG003"))
converter.convert_all()