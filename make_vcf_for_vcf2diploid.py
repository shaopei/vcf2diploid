# usuage: zcat input.vcf4.2.gz |python make_vcf_for_vcf2diploid.py output_file_name.vcf #_of_lines_in_header Col_index_of_Father Col_index_of_Mother
# example : zcat mgp.v5.merged.snps_all.dbSNP142.vcf.gz |python make_vcf_for_vcf2diploid.py CAST_129S1_F1hybrid.forAlleleSeq.vcf 104 25 11
# modified vcf4.2 to vcf4.0 file for vcf2diploid
# keep two mouse strains of interest (column) 25 11 in the example, index start at 1
# keep PASS filtered (row)
# filter out genotypes such as 0|., .|.,1|.,... only keep genotypes 1|1, 0|0, 0|1, 1|0
# filter out snps where both pareants are identical to ref genome. ie 0|0, 0|0 (child also 0|0)
# only keep snps where both parents have homolog alleles
# add artificial F1 Genotype



from sys import argv
from sys import stdin
output_f = argv[1]#'CAST_129S1_F1hybrid.forAlleleSeq.vcf'
header_line_number = int(argv[2]) #104
pat_strain = int(argv[3]) -1 #25-1 #index start at 1


print argv
chromosome_Set = set([str(i) for i in range(1,20)]+['X','Y'])

def from_two_parents():
    with open(output_f, 'w') as out:
        line_number = 1
        for l in stdin:
            if line_number > header_line_number:
                ll = l.strip().split('\t')
                if ll[6] == "PASS" and ll[0] in chromosome_Set:
                    FORMAT= ':'.join(ll[8].split(':')[0:3])
                    temp1 = ll[pat_strain].split(':')[0:3]
                    temp2 = ll[mat_strain].split(':')[0:3]
                    parent1=':'.join(['|'.join(temp1[0].split('/'))]+temp1[1:])
                    parent2=':'.join(['|'.join(temp2[0].split('/'))]+temp2[1:])
                    if parent1[0] != '.' and parent1[2] != '.' and parent2[0] != '.' and parent2[2]!= '.':
                        if parent1[0:3] != '0|0' or parent2[0:3] != '0|0' :
                            if parent1[0] == parent1[2] and parent2[0] == parent2[2]:
                                child = ':'.join([parent1[0]+'|'+parent2[0], '.','.'])
                                out.write('\t'.join(ll[0:8]+[FORMAT,parent1,parent2,child]))
                                out.write('\n')
                    #else:
                    #    print parent1,parent2
            elif line_number == header_line_number:
                ll = l.strip().split('\t')
                pat_name = '.'.join(ll[pat_strain].split('_'))
                mat_name = '.'.join(ll[mat_strain].split('_'))
                out.write('\t'.join(ll[0:9]+[pat_name,mat_name,'P.'+pat_name+'_M.'+mat_name]))
                out.write('\n')
            else:
                out.write(l)
            line_number+=1

def mother_B6():
    with open(output_f, 'w') as out:
        line_number = 1
        for l in stdin:
            if line_number > header_line_number:
                ll = l.strip().split('\t')
                if ll[6] == "PASS" and ll[0] in chromosome_Set:
                    FORMAT= ':'.join(ll[8].split(':')[0:3])
                    temp1 = ll[pat_strain].split(':')[0:3]
                    parent1=':'.join(['|'.join(temp1[0].split('/'))]+temp1[1:])
                    parent2='0|0:.:.' 
                    if parent1[0] != '.':
                        if parent1[0:3] != '0|0':
                            if parent1[0] == parent1[2]:
                                child = ':'.join([parent1[0]+'|'+parent2[0], '.','.'])
                                out.write('\t'.join(ll[0:8]+[FORMAT,parent1,parent2,child]))
                                out.write('\n')
                    #else:
                    #    print parent1
            elif line_number == header_line_number:
                ll = l.strip().split('\t')
                pat_name = '.'.join(ll[pat_strain].split('_'))
                out.write('\t'.join(ll[0:9]+[pat_name,'C57BL.6J','P.'+pat_name+'_M.C57BL.6J']))
                out.write('\n')
            else:
                out.write(l)
            line_number+=1

if __name__ == "__main__":
    if len(argv)==5:
        mat_strain = int(argv[4]) -1 #11-1 #index start at 1
        from_two_parents()
    elif len(argv)==4:
        mother_B6()
    
    
    
    