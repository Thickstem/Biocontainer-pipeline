import os
import argparse
from tqdm import tqdm
import numpy as np
import pandas as pd
from Bio import SeqIO


def _argparse():
    args = argparse.ArgumentParser()
    args.add_argument("--fasta", type=str, help="path to fasta")
    args.add_argument("--sam", type=str, help="path to input sam file")
    args.add_argument("--out_dir", type=str, help="path to dir for output file")
    opt = args.parse_args()
    return opt


def parse_fasta(fasta_path: str) -> dict:
    cds_reg_dict = {}  # {enst_id:[start_cds,end_cds]}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        enst_id = rec.id
        desc = rec.description
        seq_len = len(rec.seq)
        cat = desc.split()[2].split("=")[-1]

        if cat == "mRNA":
            cds: str = desc.split()[-1].split("=")[-1]  # cds:str = 'xx-yy'
            cds_reg_dict[enst_id] = {
                "region": [int(cds.split("-")[0]), int(cds.split("-")[1])],
                "utr5_len": int(cds.split("-")[0]),
                "cds_len": int(cds.split("-")[1]) - int(cds.split("-")[0]),
                "utr3_len": seq_len - int(cds.split("-")[1]),
            }
        else:
            cds_reg_dict[enst_id] = {
                "region": [-1, -1],
                "utr5_len": -1,
                "cds_len": -1,
                "utr3_len": -1,
            }
    return cds_reg_dict


def calc_TPM(
    df: pd.DataFrame,
    total_read_count: int,
    read_num_col: str,
    len_col: str,
    out_col: str,
):
    total_count = total_read_count
    len_per_1000bp = df[read_num_col] / df[len_col] * 10**3
    read_per_100M = len_per_1000bp / total_count * 10**6
    df[out_col] = read_per_100M

    return df


def count_dict(count_dic: dict, enst_id: str) -> dict:
    count = count_dic.get(enst_id)
    if count == None:
        count_dic[enst_id] = 1
    else:
        count_dic[enst_id] += 1

    return count_dic


def save_count_dic(
    opt: argparse.Namespace,
    fasta_df: pd.DataFrame,
    total_read_count: int,
    count_dic: dict,
    reg_name: str,
) -> None:
    count_df = pd.DataFrame(
        list(count_dic.items()), columns=["transcript_id", "read_count"]
    )
    count_df = pd.merge(count_df, fasta_df, on="transcript_id", how="inner")
    tpm_df = calc_TPM(
        count_df,
        total_read_count,
        read_num_col="read_count",
        len_col=reg_name + "_len",
        out_col="TPM_" + reg_name,
    )
    tpm_df.to_csv(os.path.join(opt.out_dir, f"ReadCount_{reg_name}.csv"))


def main(opt: argparse.Namespace):
    utr5_count_dic = {}
    cds_count_dic = {}
    utr3_count_dic = {}
    cds_reg_dict: dict = parse_fasta(opt.fasta)
    total_read_count = 0
    input_sam = open(opt.sam, "r")
    for l in tqdm(input_sam):
        if l[0] == "@":  # copy headers
            continue
        else:
            total_read_count += 1
            info = l.split()  # 0:read_id, 1:FLAG, 2:RNAME, 3:Start_pos, 4:MAPQ...
            enst_id = info[2]
            cds_reg: list = cds_reg_dict[enst_id]["region"]
            pos = int(info[3])
            ### diverge process along with mapped pos.
            if pos < cds_reg[0]:  # mapped to 5utr
                utr5_count_dic = count_dict(utr5_count_dic, enst_id)
            elif (cds_reg[0] <= pos) and (pos < cds_reg[1]):  # mapped to CDS
                cds_count_dic = count_dict(cds_count_dic, enst_id)
            elif cds_reg[1] <= pos:  # mapped to 3utr
                utr3_count_dic = count_dict(utr3_count_dic, enst_id)
            elif cds_reg[0] == -1:  # not mRNA
                continue

    fasta_df = pd.DataFrame(cds_reg_dict).T
    fasta_df = fasta_df.reset_index().rename(columns={"index": "transcript_id"})

    save_dic = {
        "utr5": utr5_count_dic,
        "cds": cds_count_dic,
        "utr3": utr3_count_dic,
    }
    for reg_name, count_dic in save_dic.items():
        save_count_dic(opt, fasta_df, total_read_count, count_dic, reg_name)


if __name__ == "__main__":
    opt = _argparse()
    main(opt)
