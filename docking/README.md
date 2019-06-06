## Flexible Docking

Flexible docking is performed using `smina` (Feb 12 2019) on the PDBbind18 dataset.

Flexible docking is performed using automatic flexible residue (`--flexdist_ligand`) and box size (`--autobox_ligand`) detection, with the following parameters (see `variables/docking`):

| Option            | Value | 
| :----------------:|:-----:|
| `--flexdist`      | 3     |
| `--autobox_add`   | 10    |
| `--exhaustiveness`| 8     |
| `--num_modes`     | 20    |

### Serial Docking

The script `docking.sh` performs docking in serial, i.e. going trough one system at a time. Given the large size of the PDBbind18 dataset and the increased computational cost, docking should be run in parallel.


### Parallel Docking

The script `docking.sub` allows to submit different docking jobs on a HPC cluster. The PDBbind18 dataset is split in chunks of 100 systems and each chunk is docked in parallel.

Only a job for the Sun Grid Engine (SGE) scheduler is provided (see `templates/sge.job`), but job scripts for other schedulers can be easily added.

### Analysis

The script `analysis.py` allows to automatically check how many systems have been docked successfully and provides some statistics. In addition, problematic systems where no conformations were found in the given search space are detected and listed.

The following `refined` systems terminate with `WARNING: Could not find any conformations completely within the search space`:
```
1k22 2bmk 3bpc 2zgx 3qtv 2jh6 1mu6 1mue 4o97 2zfp 1ype 2jh0 3f68 4hfp 4ax9
3qto 3tu7 3shc 1c4u 1nt1 1g32 1mu8 4loy 3hzm 1qbv 1g30 4e7r 1d6w 2cf8 1nm6
3si4 4isi 3qx5 2cf9 1ta6 3bv9 3hzk 3qwc 1bhx 3p17 4nj9 1zgi 4m7j 2jh5 1sl3
1c1v 1uvt 3hzv 1d9i 3sv2 1k21 2zda 1tom 1ypg 1ghy 3sha 1ypj 1z71 3utu 1sb1
2cn0 3si3 2zc9 1yei 1ghw
```

The following `other` systems terminate with `WARNING: Could not find any conformations completely within the search space`:
```
2fx8 3s7f 1qhr 1qj6 1ett 3dux 5icz 2fes 1ca8 4hs8 1wun 4zxy 1dwc 1wss 4jyu
4rn6 5wqd 2aei 2anm 1w2k 1f3j 1nzq 2bxt 3c1k 5j6d 1klg 2zp0 3chd 1o0d 1a2c
3dhk 1bmm 3tpu 4ufg 1ny2 2c8y 2ank 1wqv 1dwd 1qj7 1wtg 1c4y 1uvu 1w7x 2c8w
1awf 2r1x 2fx9 4ngh 4ufe 2l7u 2bxu 4pri 6b8i 1zgv 3egk 1uvs 1wv7 2c93 2zo3
2zzu 2a2x 1ets 2kpl 2bvx 4nga 3e0p 1sje 1bb0 1ta2 2gde 1hdt 1bmn 3v30 1c4v
1ba8 4lxb 1doj 1awh 7kme 1j81 4ufd 2c8x
```