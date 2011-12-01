import numpy

# TODO: setup script to generate config file

# Input data

project_id = 'BCSL'
project_date = '2011'
prot_db_name = 'Bacteria_prot'

genomes = [
# Finished chromosomes
{'name': 'BMB171', 'input': 'cgbk', 'file': 'CP001903.gbk',
 'acc': 'CP001903', 'size': 5330, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus thuringiensis BMB171 chromosome'},
{'name': 'CT-43', 'input': 'cgbk', 'file': 'CP001907.gbk',
 'acc': 'CP001907', 'size': 5487, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus thuringiensis sv chinensis CT-43 chromosome'},
{'name': 'YBT-020', 'input': 'cgbk', 'file': 'CP002508.gbk',
 'acc': 'CP002508', 'size': 5355, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus thuringiensis sv finitimus YBT-020 chromosome'},
{'name': 'ATCC10987', 'input': 'cgbk', 'file': 'NC_003909.gbk',
 'acc': 'NC_003909', 'size': 5224, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus cereus ATCC 10987 chromosome'},
{'name': 'ATCC14579', 'input': 'cgbk', 'file': 'NC_004722.gbk',
 'acc': 'NC_004722', 'size': 5412, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus cereus ATCC 14579 chromosome'},
{'name': '97-27', 'input': 'cgbk', 'file': 'NC_005957.gbk',
 'acc': 'NC_005957', 'size': 5238, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus thuringiensis sv konkukian str. 97-27 chromosome'},
{'name': 'E33L', 'input': 'cgbk', 'file': 'NC_006274.gbk',
 'acc': 'NC_006274', 'size': 5301, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus cereus E33L chromosome'},
{'name': 'AmesAnc', 'input': 'cgbk', 'file': 'NC_007530.gbk',
 'acc': 'NC_007530', 'size': 5227, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus anthracis str. Ames Ancestor chromosome'},
{'name': 'AlHakam', 'input': 'cgbk', 'file': 'NC_008600.gbk',
 'acc': 'NC_008600', 'size': 5257, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus thuringiensis str. Al Hakam chromosome'},
{'name': 'NVH391-98', 'input': 'cgbk', 'file': 'NC_009674.gbk',
 'acc': 'NC_009674', 'size': 4087, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus cytotoxicus NVH 391-98 chromosome'},
{'name': 'KBAB4', 'input': 'cgbk', 'file': 'NC_010184.gbk',
 'acc': 'NC_010184', 'size': 5263, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus weihenstephanensis KBAB4 chromosome'},
{'name': 'AH187', 'input': 'cgbk', 'file': 'NC_011658.gbk',
 'acc': 'NC_011658', 'size': 5269, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus cereus AH187 chromosome'},
{'name': 'B4264', 'input': 'cgbk', 'file': 'NC_011725.gbk',
 'acc': 'NC_011725', 'size': 5419, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus cereus B4264 chromosome'},
{'name': 'G9842', 'input': 'cgbk', 'file': 'NC_011772.gbk',
 'acc': 'NC_011772', 'size': 5387, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus cereus G9842 chromosome'},
{'name': 'AH820', 'input': 'cgbk', 'file': 'NC_011773.gbk',
 'acc': 'NC_011773', 'size': 5303, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus cereus AH820 chromosome'},
{'name': 'Q1', 'input': 'cgbk', 'file': 'NC_011969.gbk',
 'acc': 'NC_011969', 'size': 5214, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus cereus Q1 chromosome'},
{'name': '03BB102', 'input': 'cgbk', 'file': 'NC_012472.gbk',
 'acc': 'NC_012472', 'size': 5270, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus cereus 03BB102 chromosome'},
{'name': '03BB102', 'input': 'cgbk', 'file': 'NC_014335.gbk',
 'acc': 'NC_014335', 'size': 5196, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus cereus biovar anthracis str. CI chromosome'},
# Known pXO1-like plasmids
{'name': 'pXO1', 'input': 'cgbk', 'file': 'NC_007322.gbk',
 'acc': 'NC_007322', 'size': 182, 'offset': (0,0), 'ignore': (0, 0),
 'defline': 'Bacillus anthracis str. Ames Ancestor plasmid pXO1'},
#{'name': 'pCI-XO1', 'input': 'cgbk', 'file': 'NC_014331.gbk',
# 'acc': 'NC_014331', 'size': 182, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus biovar anthracis str. CI plasmid pCI-XO1'},
#{'name': 'pBCXO1', 'input': 'cgbk', 'file': 'NC_010934.gbk',
# 'acc': 'NC_010934', 'size': 191, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus plasmid pBCXO1'},
#{'name': 'pBc10987', 'input': 'mfas', 'file': 'AE017195.fas',
# 'acc': 'AE017195', 'size': 208, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus ATCC 10987 plasmid pBc10987'},
#{'name': 'pBc239', 'input': 'mfas', 'file': 'CP000228.fas',
# 'acc': 'CP000228', 'size': 239, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus Q1 plasmid pBc239'},
#{'name': 'pH308197_258', 'input': 'mfas', 'file': 'CP001166.fas',
# 'acc': 'CP001166', 'size': 258, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus H3081.97 plasmid pH308197_258'},
#{'name': 'pAH187_270', 'input': 'cgbk', 'file': 'NC_011655.gbk',
# 'acc': 'NC_011655', 'size': 270, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus AH187 plasmid pAH187_270'},
#{'name': 'pAH820_272', 'input': 'cgbk', 'file': 'NC_011777.gbk',
# 'acc': 'NC_011777', 'size': 272, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus AH820 plasmid pAH820_272'},
#{'name': 'p03BB102_179', 'input': 'cgbk', 'file': 'NC_012473.gbk',
# 'acc': 'NC_012473', 'size': 180, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus 03BB102 plasmid p03BB102_179'},
## Known pXO2-like plasmids
#{'name': 'pXO2', 'input': 'cgbk', 'file': 'NC_007323.gbk',
# 'acc': 'NC_007323', 'size': 95, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus anthracis str. Ames Ancestor plasmid pXO2'},
#{'name': 'pCI-XO2', 'input': 'cgbk', 'file': 'NC_014332.gbk',
# 'acc': 'NC_014332', 'size': 94, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus biovar anthracis str. CI plasmid pCI-XO2'},
#{'name': 'pAW63', 'input': 'cgbk', 'file': 'DQ025752.gbk',
# 'acc': 'DQ025752', 'size': 72, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv kurstaki str. HD73 plasmid pAW63'},
#{'name': 'pBT9727', 'input': 'cgbk', 'file': 'NC_006578.gbk',
# 'acc': 'NC_006578', 'size': 77, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv konkukian str. 97-27 plasmid pBT9727'},
## Other plasmids > 100 kb
#{'name': 'pBMB171', 'input': 'cgbk', 'file': 'CP001904.gbk',
# 'acc': 'CP001904', 'size': 313, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis BMB171 plasmid pBMB171'},
#{'name': 'pCT127', 'input': 'cgbk', 'file': 'CP001908.gbk',
# 'acc': 'CP001908', 'size': 128, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv chinensis CT-43 plasmid pCT127'},
#{'name': 'pCT281', 'input': 'cgbk', 'file': 'CP001910.gbk',
# 'acc': 'CP001910', 'size': 281, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv chinensis CT-43 plasmid pCT281'},
#{'name': 'pBMB26', 'input': 'cgbk', 'file': 'CP002509.gbk',
# 'acc': 'CP002509', 'size': 188, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv finitimus YBT-020 plasmid pBMB26'},
#{'name': 'pBMB28', 'input': 'cgbk', 'file': 'CP002510.gbk',
# 'acc': 'CP002510', 'size': 139, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv finitimus YBT-020 plasmid pBMB28'},
#{'name': 'pE33L466', 'input': 'cgbk', 'file': 'NC_007103.gbk',
# 'acc': 'NC_007103', 'size': 466, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus E33L plasmid pE33L466'},
#{'name': 'pBWB401', 'input': 'cgbk', 'file': 'NC_010180.gbk',
# 'acc': 'NC_010180', 'size': 417, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus weihenstephanensis KBAB4 plasmid pBWB401'},
#{'name': 'pG9842_140', 'input': 'cgbk', 'file': 'NC_011774.gbk',
# 'acc': 'NC_011774', 'size': 140, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus G9842 plasmid pG9842_140'},
#{'name': 'pG9842_209', 'input': 'cgbk', 'file': 'NC_011775.gbk',
# 'acc': 'NC_011775', 'size': 209, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus G9842 plasmid pG9842_209'},
## Other plasmids > 20 kb
#{'name': 'pCT51', 'input': 'cgbk', 'file': 'CP001911.gbk',
# 'acc': 'CP001911', 'size': 51, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv chinensis CT-43 plasmid pCT51'},
#{'name': 'pCT72', 'input': 'cgbk', 'file': 'CP001913.gbk',
# 'acc': 'CP001913', 'size': 72, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv chinensis CT-43 plasmid pCT72'},
#{'name': 'pCT83', 'input': 'cgbk', 'file': 'CP001915.gbk',
# 'acc': 'CP001915', 'size': 84, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv chinensis CT-43 plasmid pCT83'},
#{'name': 'pE33L54', 'input': 'cgbk', 'file': 'NC_007105.gbk',
# 'acc': 'NC_007105', 'size': 54, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus E33L plasmid pE33L54'},
#{'name': 'pALH1', 'input': 'cgbk', 'file': 'NC_008598.gbk',
# 'acc': 'NC_008598', 'size': 56, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis str. Al Hakam plasmid pALH1'},
#{'name': 'pBc53', 'input': 'cgbk', 'file': 'NC_011971.gbk',
# 'acc': 'NC_011971', 'size': 53, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus Q1 plasmid pBc53'},
#{'name': 'pBWB402', 'input': 'cgbk', 'file': 'NC_010181.gbk',
# 'acc': 'NC_010181', 'size': 75, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus weihenstephanensis KBAB4 plasmid pBWB402'},
#{'name': 'pBWB403', 'input': 'cgbk', 'file': 'NC_010182.gbk',
# 'acc': 'NC_010182', 'size': 65, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus weihenstephanensis KBAB4 plasmid pBWB403'},
#{'name': 'pBWB404', 'input': 'cgbk', 'file': 'NC_010183.gbk',
# 'acc': 'NC_010183', 'size': 53, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus weihenstephanensis KBAB4 plasmid pBWB404'},
#{'name': 'pAH187_45', 'input': 'cgbk', 'file': 'NC_011656.gbk',
# 'acc': 'NC_011656', 'size': 45, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus AH187 plasmid pAH187_45'},
## Other plasmids < 20 kb
#{'name': 'pCT6880', 'input': 'cgbk', 'file': 'CP001912.gbk',
# 'acc': 'CP001912', 'size': 7, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv chinensis CT-43 plasmid pCT6880'},
#{'name': 'pCT8252', 'input': 'cgbk', 'file': 'CP001914.gbk',
# 'acc': 'CP001914', 'size': 8, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv chinensis CT-43 plasmid pCT8252'},
#{'name': 'pCT8513', 'input': 'cgbk', 'file': 'CP001916.gbk',
# 'acc': 'CP001916', 'size': 9, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv chinensis CT-43 plasmid pCT8513'},
#{'name': 'pCT9547', 'input': 'cgbk', 'file': 'CP001917.gbk',
# 'acc': 'CP001917', 'size': 10, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv chinensis CT-43 plasmid pCT9547'},
#{'name': 'pCT14', 'input': 'cgbk', 'file': 'CP001909.gbk',
# 'acc': 'CP001909', 'size': 15, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv chinensis CT-43 plasmid pCT14'},
#{'name': 'pBC9801', 'input': 'cgbk', 'file': 'NC_009673.gbk',
# 'acc': 'NC_009673', 'size': 7, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus subsp. cytotoxis NVH 391-98 plasmid pBC9801'},
#{'name': 'pAH820_3', 'input': 'cgbk', 'file': 'NC_011776.gbk',
# 'acc': 'NC_011776', 'size': 3, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus AH820 plasmid pAH820_3'},
#{'name': 'pAH820_10', 'input': 'cgbk', 'file': 'NC_011771.gbk',
# 'acc': 'NC_011771', 'size': 11, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus AH820 plasmid pAH820_10'},
#{'name': 'pE33L8', 'input': 'cgbk', 'file': 'NC_007106.gbk',
# 'acc': 'NC_007106', 'size': 8, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus E33L plasmid pE33L8'},
#{'name': 'pE33L9', 'input': 'cgbk', 'file': 'NC_007107.gbk',
# 'acc': 'NC_007107', 'size': 9, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus E33L plasmid pE33L9'},
#{'name': 'pE33L5', 'input': 'cgbk', 'file': 'NC_007104.gbk',
# 'acc': 'NC_007104', 'size': 5, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus E33L plasmid pE33L5'},
#{'name': 'pAH187_3', 'input': 'cgbk', 'file': 'NC_011657.gbk',
# 'acc': 'NC_011657', 'size': 3, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus AH187 plasmid pAH187_3'},
#{'name': 'pAH187_12', 'input': 'cgbk', 'file': 'NC_011654.gbk',
# 'acc': 'NC_011654', 'size': 12, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus AH187 plasmid pAH187_12'},
#{'name': 'pBAslCI14', 'input': 'cgbk', 'file': 'NC_014333.gbk',
# 'acc': 'NC_014333', 'size': 14, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus biovar anthracis str. CI plasmid pBAslCI14'},
#{'name': 'pBClin15', 'input': 'cgbk', 'file': 'NC_004721.gbk',
# 'acc': 'NC_004721', 'size': 15, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus ATCC 14579 plasmid pBClin15'},
## Genbank WGS genomes
#{'name': 'G9241', 'input': 'mfas', 'file': 'NZ_AAEK0.fas',
# 'acc': 'NZ_AAEK0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus G9241'},
#{'name': 'ATCC35646', 'input': 'mfas', 'file': 'NZ_AAJM0.fas',
# 'acc': 'NZ_AAJM0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis serovar israelensis ATCC 35646'},
#{'name': 'AH1134', 'input': 'mfas', 'file': 'NZ_ABDA0.fas',
# 'acc': 'NZ_ABDA0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus AH1134'},
#{'name': 'H308197', 'input': 'mfas', 'file': 'NZ_ABDL0.fas',
# 'acc': 'NZ_ABDL0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus H3081.97'},
#{'name': '03BB108', 'input': 'mfas', 'file': 'NZ_ABDM0.fas',
# 'acc': 'NZ_ABDM0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus 03BB108'},
#{'name': 'm1293', 'input': 'mfas', 'file': 'NZ_ACLS0.fas',
# 'acc': 'NZ_ACLS0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus m1293'},
#{'name': 'ATCC10876', 'input': 'mfas', 'file': 'NZ_ACLT0.fas',
# 'acc': 'NZ_ACLT0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus ATCC 10876'},
#{'name': 'BGSC6E1', 'input': 'mfas', 'file': 'NZ_ACLU0.fas',
# 'acc': 'NZ_ACLU0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BGSC 6E1'},
#{'name': '172560W', 'input': 'mfas', 'file': 'NZ_ACLV0.fas',
# 'acc': 'NZ_ACLV0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus 172560W'},
#{'name': 'MM3', 'input': 'mfas', 'file': 'NZ_ACLW0.fas',
# 'acc': 'NZ_ACLW0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus MM3'},
#{'name': 'AH621', 'input': 'mfas', 'file': 'NZ_ACLX0.fas',
# 'acc': 'NZ_ACLX0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus AH621'},
#{'name': 'R309803', 'input': 'mfas', 'file': 'NZ_ACLY0.fas',
# 'acc': 'NZ_ACLY0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus R309803'},
#{'name': 'ATCC4342', 'input': 'mfas', 'file': 'NZ_ACLZ0.fas',
# 'acc': 'NZ_ACLZ0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus ATCC 4342'},
#{'name': 'm1550', 'input': 'mfas', 'file': 'NZ_ACMA0.fas',
# 'acc': 'NZ_ACMA0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus m1550'},
#{'name': 'BDRD-ST24', 'input': 'mfas', 'file': 'NZ_ACMB0.fas',
# 'acc': 'NZ_ACMB0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BDRD-ST24'},
#{'name': 'BDRD-ST26', 'input': 'mfas', 'file': 'NZ_ACMC0.fas',
# 'acc': 'NZ_ACMC0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BDRD-ST26'},
#{'name': 'BDRD-ST196', 'input': 'mfas', 'file': 'NZ_ACMD0.fas',
# 'acc': 'NZ_ACMD0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BDRD-ST196'},
#{'name': 'BDRD-Cer4', 'input': 'mfas', 'file': 'NZ_ACME0.fas',
# 'acc': 'NZ_ACME0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BDRD-Cer4'},
#{'name': '95/8201', 'input': 'mfas', 'file': 'NZ_ACMF0.fas',
# 'acc': 'NZ_ACMF0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus 95/8201'},
#{'name': 'Rock1-3', 'input': 'mfas', 'file': 'NZ_ACMG0.fas',
# 'acc': 'NZ_ACMG0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus Rock1-3'},
#{'name': 'Rock1-15', 'input': 'mfas', 'file': 'NZ_ACMH0.fas',
# 'acc': 'NZ_ACMH0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus Rock1-15'},
#{'name': 'Rock3-28', 'input': 'mfas', 'file': 'NZ_ACMI0.fas',
# 'acc': 'NZ_ACMI0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus Rock3-28'},
#{'name': 'Rock3-29', 'input': 'mfas', 'file': 'NZ_ACMJ0.fas',
# 'acc': 'NZ_ACMJ0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus Rock3-29'},
#{'name': 'Rock3-42', 'input': 'mfas', 'file': 'NZ_ACMK0.fas',
# 'acc': 'NZ_ACMK0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus Rock3-42'},
#{'name': 'Rock3-44', 'input': 'mfas', 'file': 'NZ_ACML0.fas',
# 'acc': 'NZ_ACML0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus Rock3-44'},
#{'name': 'Rock4-2', 'input': 'mfas', 'file': 'NZ_ACMM0.fas',
# 'acc': 'NZ_ACMM0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus Rock4-2'},
#{'name': 'Rock4-18', 'input': 'mfas', 'file': 'NZ_ACMN0.fas',
# 'acc': 'NZ_ACMN0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus Rock4-18'},
#{'name': 'F65185', 'input': 'mfas', 'file': 'NZ_ACMO0.fas',
# 'acc': 'NZ_ACMO0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus F65185'},
#{'name': 'AH603', 'input': 'mfas', 'file': 'NZ_ACMP0.fas',
# 'acc': 'NZ_ACMP0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus AH603'},
#{'name': 'AH676', 'input': 'mfas', 'file': 'NZ_ACMQ0.fas',
# 'acc': 'NZ_ACMQ0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus AH676'},
#{'name': 'AH1271', 'input': 'mfas', 'file': 'NZ_ACMR0.fas',
# 'acc': 'NZ_ACMR0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus AH1271'},
#{'name': 'AH1272', 'input': 'mfas', 'file': 'NZ_ACMS0.fas',
# 'acc': 'NZ_ACMS0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus AH1272'},
#{'name': 'AH1273', 'input': 'mfas', 'file': 'NZ_ACMT0.fas',
# 'acc': 'NZ_ACMT0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus AH1273'},
#{'name': 'DSM 2048', 'input': 'mfas', 'file': 'NZ_ACMU0.fas',
# 'acc': 'NZ_ACMU0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus mycoides DSM 2048'},
#{'name': 'Rock1-4', 'input': 'mfas', 'file': 'NZ_ACMV0.fas',
# 'acc': 'NZ_ACMV0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus mycoides Rock1-4'},
#{'name': 'Rock3-17', 'input': 'mfas', 'file': 'NZ_ACMW0.fas',
# 'acc': 'NZ_ACMW0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus mycoides Rock3-17'},
#{'name': 'DSM12442', 'input': 'mfas', 'file': 'NZ_ACMX0.fas',
# 'acc': 'NZ_ACMX0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus pseudomycoides DSM 12442'},
#{'name': 'BGSC4Y1', 'input': 'mfas', 'file': 'NZ_ACMY0.fas',
# 'acc': 'NZ_ACMY0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BGSC 4Y1'},
#{'name': 'Bt407', 'input': 'mfas', 'file': 'NZ_ACMZ0.fas',
# 'acc': 'NZ_ACMZ0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bt407'},
#{'name': 'T01001', 'input': 'mfas', 'file': 'NZ_ACNA0.fas',
# 'acc': 'NZ_ACNA0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv thuringiensis str. T01001'},
#{'name': 'T04001', 'input': 'mfas', 'file': 'NZ_ACNB0.fas',
# 'acc': 'NZ_ACNB0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv sotto str. T04001'},
#{'name': 'T13001', 'input': 'mfas', 'file': 'NZ_ACNC0.fas',
# 'acc': 'NZ_ACNC0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv pakistani str. T13001'},
#{'name': 'T03a001', 'input': 'mfas', 'file': 'NZ_ACND0.fas',
# 'acc': 'NZ_ACND0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis sv kurstaki str. T03a001'},
#{'name': 'BGSC4AJ1', 'input': 'mfas', 'file': 'NZ_ACNE0.fas',
# 'acc': 'NZ_ACNE0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BGSC 4AJ1'},
#{'name': 'ATCC10792', 'input': 'mfas', 'file': 'NZ_ACNF0.fas',
# 'acc': 'NZ_ACNF0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'berliner ATCC 10792'},
#{'name': 'BGSC4AW1', 'input': 'mfas', 'file': 'NZ_ACNG0.fas',
# 'acc': 'NZ_ACNG0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BGSC 4AW1'},
#{'name': 'BGSC4BA1', 'input': 'mfas', 'file': 'NZ_ACNH0.fas',
# 'acc': 'NZ_ACNH0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BGSC 4BA1'},
#{'name': 'BGSC4BD1', 'input': 'mfas', 'file': 'NZ_ACNI0.fas',
# 'acc': 'NZ_ACNI0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BGSC 4BD1'},
#{'name': 'BGSC4CC1', 'input': 'mfas', 'file': 'NZ_ACNJ0.fas',
# 'acc': 'NZ_ACNJ0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BGSC 4CC1'},
#{'name': 'IBL200', 'input': 'mfas', 'file': 'NZ_ACNK0.fas',
# 'acc': 'NZ_ACNK0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis IBL 200'},
#{'name': 'IBL4222', 'input': 'mfas', 'file': 'NZ_ACNL0.fas',
# 'acc': 'NZ_ACNL0', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus thuringiensis IBL 4222'},
#{'name': 'SJ1', 'input': 'mfas', 'file': 'NZ_ADFM0.fas',
# 'acc': 'NZ_ADFM0',  'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus SJ1'},
## Broad genomes
#{'name': 'IS075', 'input': 'mfas', 'file': 'IS075.fas',
# 'acc': 'IS075', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus IS075 (Broad batch 1)'},
#{'name': 'Schrouff', 'input': 'mfas', 'file': 'Schrouff.fas',
# 'acc': 'Schrouff', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus Schrouff contigs (Broad batch 1)'},
#{'name': 'TIAC129', 'input': 'mfas', 'file': 'TIAC129.fas',
# 'acc': 'TIAC129', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus TIAC129 contigs (Broad batch 1)'},
#{'name': 'VD142', 'input': 'mfas', 'file': 'VD142.fas',
# 'acc': 'VD142', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus VD142 contigs (Broad batch 1)'},
#{'name': 'VD022', 'input': 'mfas', 'file': 'VD022.fas',
# 'acc': 'VD022', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus VD022 contigs (Broad batch 1)'},
#{'name': 'VD022_454', 'input': 'mfas', 'file': 'VD022_454.fas',
# 'acc': 'VD022_454', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'Bacillus cereus VD022 contigs (454 version)'},
#{'name': 'BAG1X1-1a', 'input': 'mfas', 'file': 'G13150_assembly.fasta',
# 'acc': 'BAG1X1-1a', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BAG1X1-1 assembly  (Broad batch 2)'},
#{'name': 'BAG1X1-1c', 'input': 'mfas', 'file': 'G13150_contigs.fasta',
# 'acc': 'BAG1X1-1c', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BAG1X1-1 contigs  (Broad batch 2)'},
#{'name': 'BAG1X1-2a', 'input': 'mfas', 'file': 'G13151_assembly.fasta',
# 'acc': 'BAG1X1-2a', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BAG1X1-2 assembly  (Broad batch 2)'},
#{'name': 'BAG1X1-2c', 'input': 'mfas', 'file': 'G13151_contigs.fasta',
# 'acc': 'BAG1X1-2c', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BAG1X1-2 contigs  (Broad batch 2)'},
#{'name': 'BAG1X2-1a', 'input': 'mfas', 'file': 'G13156_assembly.fasta',
# 'acc': 'BAG1X2-1a', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BAG1X2-1 assembly  (Broad batch 2)'},
#{'name': 'BAG1X2-1c', 'input': 'mfas', 'file': 'G13156_contigs.fasta',
# 'acc': 'BAG1X2-1c', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BAG1X2-1 contigs  (Broad batch 2)'},
#{'name': 'BAG1O-2a', 'input': 'mfas', 'file': 'G13159_assembly.fasta',
# 'acc': 'BAG1O-2a', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BAG1O-2 assembly  (Broad batch 2)'},
#{'name': 'BAG1O-2c', 'input': 'mfas', 'file': 'G13159_contigs.fasta',
# 'acc': 'BAG1O-2c', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BAG1O-2 contigs  (Broad batch 2)'},
#{'name': 'BAG1O-3a', 'input': 'mfas', 'file': 'G13160_assembly.fasta',
# 'acc': 'BAG1O-3a', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BAG1O-3 assembly  (Broad batch 2)'},
#{'name': 'BAG1O-3c', 'input': 'mfas', 'file': 'G13160_contigs.fasta',
# 'acc': 'BAG1O-3c', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BAG1O-3 contigs  (Broad batch 2)'},
#{'name': 'BAG1X1-3a', 'input': 'mfas', 'file': 'G13172_assembly.fasta',
# 'acc': 'BAG1X1-3a', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BAG1X1-3 assembly  (Broad batch 2)'},
#{'name': 'BAG1X1-3c', 'input': 'mfas', 'file': 'G13172_contigs.fasta',
# 'acc': 'BAG1X1-3c', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BAG1X1-3 contigs  (Broad batch 2)'},
#{'name': 'AND1407a', 'input': 'mfas', 'file': 'G13175_assembly.fasta',
# 'acc': 'AND1407a', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'AND1407 assembly  (Broad batch 2)'},
#{'name': 'AND1407c', 'input': 'mfas', 'file': 'G13175_contigs.fasta',
# 'acc': 'AND1407c', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'AND1407 contigs  (Broad batch 2)'},
#{'name': 'BAG2X1-3a', 'input': 'mfas', 'file': 'G13211_assembly.fasta',
# 'acc': 'BAG2X1-3a', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BAG2X1-3 assembly  (Broad batch 2)'},
#{'name': 'BAG2X1-3c', 'input': 'mfas', 'file': 'G13211_contigs.fasta',
# 'acc': 'BAG2X1-3c', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'BAG2X1-3 contigs  (Broad batch 2)'},
#{'name': 'VD_022a', 'input': 'mfas', 'file': 'G8177_assembly.fasta',
# 'acc': 'VD_022a', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'VD_022 assembly  (Broad batch 2)'},
#{'name': 'VD_022c', 'input': 'mfas', 'file': 'G8177_contigs.fasta',
# 'acc': 'VD_022c', 'size': 0, 'offset': (0,0), 'ignore': (0, 0),
# 'defline': 'VD_022 contigs  (Broad batch 2)'}
]

references = [#{'name': 'pXO1', 'file': 'pXO1.gbk', 'input': 'gbk',
#               'seg_mode': 'list', 'segs':
#               ({'coords': (1127,6561), 'name': 'A', 'note': '3_10'},
#                {'coords': (8991,13756),'name': 'B', 'note': '12_17'},
#                {'coords': (19200,20894),'name': 'C', 'note': '20_20'},
#                {'coords': (21714,25117),'name': 'D', 'note': '23_25'},
#                {'coords': (28490,30979),'name': 'E', 'note': '30_34'},
#                {'coords': (38284,40203),'name': 'F', 'note': '45_46'},
#                {'coords': (52063,55796),'name': 'G', 'note': '64_65'},
#                {'coords': (55987,60586),'name': 'H', 'note': '66_69'},
#                {'coords': (61404,63518),'name': 'I', 'note': '71_74'},
#                {'coords': (64090,68866),'name': 'J', 'note': '75_79'},
#                {'coords': (69021,74102),'name': 'K', 'note': '80_84'},
#                {'coords': (74119,77871),'name': 'L', 'note': '85_90'},
#                {'coords': (79735,82867),'name': 'M', 'note': '93_97'},
#                {'coords': (83606,86696),'name': 'N', 'note': '98_101'},
#                {'coords': (87297,91464),'name': 'O', 'note': '103_107'},
#                {'coords': (91495,95163),'name': 'P', 'note': '108_108'},
#                {'coords': (95212,98847),'name': 'Q', 'note': '109_114'},
#                {'coords': (100828,102495),'name': 'R', 'note': '117_118'},
#                {'coords': (155257,157413),'name': 'S', 'note': '182_184'},
#                {'coords': (174223,176584),'name': 'T', 'note': '207_212'})},
              {'name': 'pXO1', 'file': 'NC_007322.gbk', 'input': 'gbk',
               'seg_mode': 'chop', 'chop_size': 2000}]

# Function categories and legend (keys MUST be lowercase)

fct_flags = {'mge': ('transposase', 'transposon', 'intron',
                     'tyrosine recombinase', 'dna-invertase',
                     'reverse transcriptase'),
             'rep': ('something else', 'replication'),
             'syn': ('synthase', 'something else'),
             'tox': ('protective antigen', 'lethal factor',
                     'virulence factor'),
             'ger': ('spore germination protein', 'something else'),
             'tra': ('type iv', 'type ii/iv', 'topoisomerase',
                     'dna translocase ftsk', 'conjugal transfer',
                     'conjugation protein'),
             'ctl': ('transcriptional regulator', 'regulatory protein',
                     'regulator', 'transcriptional repressor'
                     'response regulator aspartate phosphatase'),
             'unk': ('uncharacterized', 'conserved domain protein'),
             'def': ('no match', 'hypothetical protein', 'pXO1-')} # default

fct_colors = {'mge': ('#66CC00', 'MGE'),
              'rep': ('#FF9900', 'Replication'),
              'syn': ('#FF00CC', 'Synthesis'),
              'tox': ('#FF0000', 'Pathogenesis'),
              'ger': ('#993333', 'Germination'),
              'tra': ('#6666FF', 'Transfer'),
              'ctl': ('#FFCC00', 'Control'),
              'unk': ('#666666', 'Uncharacterized'),
              'oth': ('#CCCCCC', 'Other'),
              'def': ('#FFFFFF', 'No match')}

# Concat contigs separator
separator = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"

# Project root directory
p_root_dir = 'data/'+project_id+'/'

# run-independent directories
fixed_dirs = {
'ref_dbs_dir': 'data/ref_dbs/',
'ori_g_dir': p_root_dir+'genomes/original/',
'gbk_contigs_dir': p_root_dir+'genomes/genbank_ctgs/',
'fas_contigs_dir': p_root_dir+'genomes/fasta_ctgs/',
'mfas_contigs_dir': p_root_dir+'genomes/mfasta_ctgs/',
'blast_db_dir': p_root_dir+'genomes/blast_db/',
'annot_trn_dir': p_root_dir+'genomes/annotation/training/',
'ctg_cds_dir': p_root_dir+'genomes/annotation/genes/',
'ctg_prot_dir': p_root_dir+'genomes/annotation/proteins/',
'ctg_blast_dir': p_root_dir+'genomes/annotation/blastp/',
'ctg_stats': p_root_dir+'genomes/contig_stats/'
}

# run-dependent directories
run_dirs = {
'ref_pickles': 'references/',
'ref_seg_dir': 'references/segments/',
'ref_gbk_dir': 'references/genbank/',
'ref_fas_dir': 'references/fasta/',
'ref_map_dir': 'maps/references/',
'match_pickles': 'matching/',
'blast_out_dir': 'matching/blastn/',
'match_out_dir': 'matching/matches/',
'scaffolds_dir': 'scaffolds/',
'run_gbk_ctgs_dir': 'gbk_contigs/',
'mauve_out_dir': 'alignments/mauve_out/',
'aln_seg_dir': 'alignments/aln_segments/',
'maps_dir': 'maps/',
'reports': 'reports/'
}

# Blast parameters
blast_prefs = {'evalue': 0.01,
               'outfmt_pref': 6}
min_match = 500     # min size for a blast hit to be considered relevant
min_score = 1000    # min score for a blast hit to be considered relevant

# Blast results arrays datatypes
blast_dtypes = numpy.dtype([('query', 'S16'),
                           ('dbhit', 'S32'),
                           ('idp', 'float'),
                           ('mlen', 'uint8'),
                           ('mms', 'uint8'),
                           ('gaps', 'uint8'),
                           ('q_start', 'uint16'),
                           ('q_end', 'uint16'),
                           ('r_start', 'uint16'),
                           ('r_end', 'uint16'),
                           ('evalue', 'S5'),
                           ('bitscore', 'float')])

# Mauve executable
mauve_exec = 'progressiveMauve'

# Datatypes
mtype = [('A','i4'), ('B','i4'), ('C', 'i4'), ('D','i4')]
segtype = [('A','i4'), ('B','i4'), ('C', 'i4'), ('D','i4'), ('E','i4')]

# Proximity thresholds for clumping
prox_D = 2000   # for ballpark estimation
prox_F = 100    # for fine alignment (set to 0 to skip clumping)

# Chopping size to limit length of detailed alignments
max_size = 3000
chop_mode = 'exact_size' # 'maxsize_bisect','maxsize_divisor','count_divisor'

# Identity percentage cutoffs and color coding
idpt = {95: '#444444',     # top similarity class (HexColor('#444444'))
        85: '#777777',     # upper middle class (HexColor('#777777'))
        70: '#BBBBBB',     # lower middle class (HexColor('#BBBBBB'))
        50: '#DDDDDD',     # low similarity class (HexColor('#DDDDDD'))
         0: '#FFFFFF'}     # lower than cutoff (HexColor('#FFFFFF'))

# Thresholds for binning contig sizes
ctg_thresholds = [100, 500, 1000]
