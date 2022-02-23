# Test data

### Edge cases to test

Case A: Default (GATC only)
Case B: Hard-masking (Ns)
Case C: Degenerate Nucleotides (=MRSVWYHKDB)
Case D: Soft-masking (gatc or mrsvwyhkdb)
Case BC: Hard-masking & Degenerate Nucleotides
Case BD: Hard-masking & Soft-masking
Case CD: Degenerate Nucleotides & Soft-masking
Case BCD: Hard-masking & Degenerate Nucleotides & Soft-masking


| File     | Base Case  | Mod A | Mod B | Mod C | Mod D | Mod BC | Mod BD | Mod CD | Mod BCD |
|----------|------------|-------|-------|-------|-------|--------|--------|--------|---------|
| test1.fa | A          | A_A   | B_A   | C_A   | D_A   | BC_A   | BD_A   | CD_A   | BCD_A   |
| test2.fa | B          | A_B   | B_B   | C_B   | D_B   | BC_B   | BD_B   | CD_B   | BCD_B   |
| test3.fa | C          | A_C   | B_C   | C_C   | D_C   | BC_C   | BD_C   | CD_C   | BCD_C   |
| test4.fa | D          | A_D   | B_D   | C_D   | D_D   | BC_D   | BD_D   | CD_D   | BCD_D   |
| test5.fa | BC         | A_BC  | B_BC  | C_BC  | D_BC  | BC_BC  | BD_BC  | CD_BC  | BCD_BC  |
| test6.fa | BD         | A_BD  | B_BD  | C_BD  | D_BD  | BC_BD  | BD_BD  | CD_BD  | BCD_BD  |
| test7.fa | CD         | A_CD  | B_CD  | C_CD  | D_CD  | BC_CD  | BD_CD  | CD_CD  | BCD_CD  |
| test8.fa | BCD        | A_BCD | B_BCD | C_BCD | D_BCD | BC_BCD | BD_BCD | CD_BCD | BCD_BCD |