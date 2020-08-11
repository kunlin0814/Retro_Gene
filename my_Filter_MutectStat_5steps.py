#!/usr/bin/python

import sys

if __name__ == '__main__':
	file0 = sys.argv[1]
	pass0 = sys.argv[2]
	whyout = sys.argv[3]

	with open(file0,'r') as f:
		file = f.read()
	with open(pass0,'r') as f:
		pas = f.read()
	lst = file.split('\n')[:-1]
	lst2 = pas.split('\n')[:-1]
	# output
	out = open(pass0 + '_filteredMut','w')
	whyout_output = open(whyout,'w')
	dic = {}
	failed = []
	for i in range(len(lst)):
		info = lst[i].split('\t')
		chrom = info[0]
		pos = info[1]
		tRef = int(info[2])
		tAlt = int(info[3])
		nRef = int(info[4])
		nAlt = int(info[5])
		total_tumor_depth = tRef + tAlt # total tumor read depth
		vaf = float(tAlt) / total_tumor_depth # VAF for tumor
		if total_tumor_depth < 10:
			failed.append(i)
			whyout_output.write(str(chrom)+'_'+str(pos)+'\t'+ "total_tumor_depth < 10 and the death is ," + str(total_tumor_depth)+"\n")
			
		else:
			if vaf < 0.05:
				failed.append(i)
				whyout_output.write(str(chrom)+'_'+str(pos)+'\t'+"VAF < 0.05 and the VAF is ,"+ str(vaf)+ '\n')

			else:
				if tAlt <= 5 and vaf < 0.15:
					failed.append(i)
					whyout_output.write(str(chrom)+'_'+str(pos)+'\t'+"tAlt <= 5 and vaf < 0.15, and tAlt and vaf is "+ str(tAlt)+','+str(vaf)+"\n")
				else:
					if total_tumor_depth < 20 and vaf < 0.2:
						failed.append(i)
						whyout_output.write(str(chrom)+'_'+str(pos)+'\t'+"total_tumor_depth < 20 and vaf < 0.2 and total_tumor_depth and vaf is "+ str(total_tumor_depth)+","+str(vaf)+"\n")
					else:
						if nAlt >= 3:
							failed.append(i)
							whyout_output.write(str(chrom)+'_'+str(pos)+'\t'+"nAlt >= 3 and nAlt is "+ str(nAlt)+'\n')
						else:
							if chrom not in dic.keys():
								dic[chrom] = []
							dic[chrom].append(pos)
	for i in range(len(lst2)):
		info = lst2[i].split('\t')
		chrom = info[0]
		pos = info[1]
		if chrom in dic.keys():
			if pos in dic[chrom]:
				string = lst2[i] + '\n'
				out.write(string)
	out.close()
