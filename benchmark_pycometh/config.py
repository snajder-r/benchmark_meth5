class ModuleConfig:
    def __init__(self):
        self.samples = ["HG002", "HG003", "HG004"]
        self.reference = "/home/r933r/data/resource/human/hg38/Ensembl_101/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        self.gff_file = "/home/r933r/data/projects/nanopore/pycometh_benchmark/annotation/annotation.gff3"
        self.methylkit_segmentation = {
            "HG003_R1": "/home/r933r/data/projects/nanopore/pycometh_benchmark/bs/HG003_rep1_seg.bed",
            "HG003_R2": "/home/r933r/data/projects/nanopore/pycometh_benchmark/bs/HG003_rep2_seg.bed",
            "HG003_mock_from_np": "/home/r933r/data/projects/nanopore/pycometh_benchmark/bs/HG003_mockbsseq_seg.bed",
            "HG003_mock_from_np_hp": "/home/r933r/data/projects/nanopore/pycometh_benchmark/bs/HG003_mockbsseq_hps_fastseg.bed",
            "HG003": "/home/r933r/data/projects/nanopore/pycometh_benchmark/bs/HG003_rep1_rep2_seg.bed",
            "HG004": "/home/r933r/data/projects/nanopore/pycometh_benchmark/bs/HG004_rep1_rep2_seg.bed",
            "parents_mock_from_np_hp": "/home/r933r/data/projects/nanopore/pycometh_benchmark/bs/parents_mockbsseq_hps_fastseg.bed",
            "parents_mock_from_np_diffmet": "/omics/groups/OE0540/internal/projects/nanopore/pycometh_benchmark/bs/parents_mockbsseq_diff_methylkit_seg.bed",
        }
        self.nanoepiseg_segmentation = {
            "HG003": "/home/r933r/data/projects/nanopore/pycometh_benchmark/segmentation/HG003.bed",
            "HG003_35": "/home/r933r/data/projects/nanopore/pycometh_benchmark/segmentation/HG003_35.bed",
            "parents": "/home/r933r/data/projects/nanopore/pycometh_benchmark/segmentation_parents/segmentation_parents.bed",
            "parents_35": "/home/r933r/data/projects/nanopore/pycometh_benchmark/segmentation_parents/segmentation_parents_mindiff_35.bed",
        }
        self.meth5_template_file = (
            "/home/r933r/data/projects/nanopore/pycometh_benchmark/met_merged/{sample}_cpg_lzfcompressed.h5"
        )
        self.meth5_gzipped_template_file = (
            "/home/r933r/data/projects/nanopore/pycometh_benchmark/met_merged/{" "sample}_cpg.h5"
        )
        self.haplotype_mapping_file = "/home/r933r/data/projects/nanopore/pycometh_benchmark/svs/haplotype_mapping.tsv"
        self.mock_bsseq_template_file = (
            "/home/r933r/data/projects/nanopore/pycometh_benchmark/met_merged/{sample}_mock_bsseq.bedGraph"
        )
        self.mock_bsseq_template_file_consistency_testing = "/home/r933r/data/projects/nanopore/pycometh_benchmark/met_merged/{sample}_mock_bsseq_small/subset{iteration}.bedGraph"
        self.mock_bsseq_template_file_consistency_methylkit_result = (
            "/home/r933r/data/projects/nanopore/pycometh_benchmark/bs/{sample}_small_example_consistency.tsv"
        )
        self.mock_bsseq_hp_template_file = (
            "/home/r933r/data/projects/nanopore/pycometh_benchmark/met_merged/{sample}_{hp}_mock_bsseq.bedGraph"
        )
        self.vcfs = {
            "HG002": "/home/r933r/data/projects/nanopore/pycometh_benchmark/svs/HG002_phased.vcf.gz",
            "HG003": "/home/r933r/data/projects/nanopore/pycometh_benchmark/svs/HG003_phased.vcf.gz",
            "HG004": "/home/r933r/data/projects/nanopore/pycometh_benchmark/svs/HG004_phased.vcf.gz",
        }
        
        self.bs_seq_files = {
            "HG003": {
                "R1": "/home/r933r/data/projects/nanopore/data/giab/bs_seq_seqc2/GSM5649434_TruSeq_HG003_LAB01_REP01.bedGraph.gz",
                "R2": "/home/r933r/data/projects/nanopore/data/giab/bs_seq_seqc2/GSM5649433_TruSeq_HG003_LAB01_REP02.bedGraph.gz",
                "combined": "/home/r933r/data/projects/nanopore/pycometh_benchmark/bs/HG003_rep1_rep2_merged.bedGraph.gz",
            },
            "HG004": {
                "R1": "/home/r933r/data/projects/nanopore/data/giab/bs_seq_seqc2/GSM5649432_TruSeq_HG004_LAB01_REP01.bedGraph.gz",
                "R2": "/home/r933r/data/projects/nanopore/data/giab/bs_seq_seqc2/GSM5649430_TruSeq_HG004_LAB01_REP02.bedGraph.gz",
                "combined": "/home/r933r/data/projects/nanopore/pycometh_benchmark/bs/HG004_rep1_rep2_merged.bedGraph.gz",
            },
        }
        self.subset_m5_file = "/home/r933r/data/projects/nanopore/pycometh_benchmark/subset_m5/{sample}_cpg.h5"
        self.modbam_template_file = (
            "/home/r933r/data/projects/nanopore/pycometh_benchmark/met_modbam_npdefault/{sample}_cpg.sorted.bam"
        )
        self.modcram_template_file = (
            "/home/r933r/data/projects/nanopore/pycometh_benchmark/met_modbam_npdefault/{" "sample}_cpg.sorted.cram"
        )
    
    def __item__(self, key):
        if key in dir(self):
            return getattr(self, key)


module_config = ModuleConfig()
