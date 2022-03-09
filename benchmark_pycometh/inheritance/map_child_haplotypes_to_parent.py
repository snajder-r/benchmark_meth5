import gzip

import pandas as pd

from benchmark_pycometh.config import module_config

def read_sv_file(path):
    with gzip.open(path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            line = line.split("\t")
            hp = line[9][:3]
            if hp[1] == "/":
                # unphased
                continue
            if hp[0] == hp[2]:
                # heterozygous
                continue
            hp = 1 if hp[0] == "1" else 2
            yield dict(chrom=line[0], pos=int(line[1]), ref=line[3], alt=line[4], hp=hp)
            

def load_svs():
    child_df = pd.DataFrame(read_sv_file(module_config.vcfs["HG002"])).set_index(["chrom", "pos", "alt"])
    pa_df = pd.DataFrame(read_sv_file(module_config.vcfs["HG003"])).set_index(["chrom", "pos", "alt"])
    pb_df = pd.DataFrame(read_sv_file(module_config.vcfs["HG004"])).set_index(["chrom", "pos", "alt"])
    return child_df, pa_df, pb_df


if __name__ == "__main__":
    child_df, pa_df, pb_df = load_svs()
    
    child_df.set_index(["chrom", "pos", "alt"]).join()
    
    
    merged_A = child_df.join(pa_df, how="inner", lsuffix="", rsuffix="_p")
    merged_B = child_df.join(pb_df, how="inner", lsuffix="", rsuffix="_p")
    
    in_both = merged_A.index.intersection(merged_B.index)
    
    merged_A = merged_A.loc[merged_A.index.difference(in_both)]
    merged_B = merged_B.loc[merged_B.index.difference(in_both)]
    
    merged_A = merged_A.reset_index()
    merged_B = merged_B.reset_index()
    
    merged_A["source"] = "HG003"
    merged_B["source"] = "HG004"
    
    merged = pd.concat([merged_A, merged_B]).sort_values(["hp", "chrom", "pos", "source"])
    merged["source_hp"] = merged.apply(lambda x: f"{x['source']}_{x['hp_p']}", axis=1)
    merged = merged.groupby("hp")
    
    
    segments = []
    for hp in [1,2]:
        start = 0
        end = 0
        source = None
        chrom = None
        for row in merged.get_group(hp).itertuples():
            if source is None:
                source = row.source_hp
                chrom = row.chrom
            
            if source != row.source_hp or chrom != row.chrom:
                if chrom == row.chrom:
                    end = row.pos
                source = source.split("_")
                segments.append(dict(chrom=chrom,start=start,end=end,child_hp=hp,parent=source[0], parent_hp=source[1]))
                start = row.pos
                chrom = row.chrom
                source = row.source_hp
            end = row.pos
    
    mapping_hps = pd.DataFrame(segments).sort_values(["chrom", "start"]).reset_index(drop=True)
    mapping_hps.to_csv(module_config.haplotype_mapping_file, sep="\t", index=False)
    
    num_strange_mappings = 0
    """ SANITY CHECKING since we don't expect the same source haplotype to provide both child haplotypes"""
    lastrow = None
    for row in mapping_hps.sort_values(["chrom", "parent", "parent_hp", "start"]).itertuples():
        if lastrow is None:
            lastrow = row
            continue
        if row.parent != lastrow.parent or row.parent_hp != lastrow.parent_hp or row.chrom != lastrow.chrom:
            lastrow = row
            continue
        
        # same source
        if row.start < lastrow.end and lastrow.start < row.end:
            if row.child_hp != lastrow.child_hp:
                num_strange_mappings+=2
                print(f"Parent {row.parent} {row.parent_hp} provided both child HPs in ", row.chrom, row.start, row.end, lastrow.start, lastrow.end)
    print(f"Number of mappings with overlaps: {num_strange_mappings} ({100*num_strange_mappings/mapping_hps.shape[0]}%)")
            