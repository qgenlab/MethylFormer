import os,sys

def write_cpg( mw, cur_chr, ful_seq, cpg_id):
   for _i in range(len(ful_seq)):
      if ful_seq[_i:_i+2]=='CG':
         mw.write("{}\t{}\t{}\t{}\n".format( cur_chr, _i,_i+1, '+' ))
         cpg_id = cpg_id + 1
   return cpg_id

if __name__=='__main__':
   ref_file =sys.argv[1]
   op_file = sys.argv[2]

   orderedchr = []
   orderedchrnum = []
   for ic in range(1, 100):
      orderedchr.append('chr'+str(ic));
      orderedchrnum.append(str(ic));

   orderedchr.append('chrX'); orderedchrnum.append('X')
   orderedchr.append('chrY'); orderedchrnum.append('Y')
   orderedchr.append('chrM'); orderedchrnum.append('M')

   ref_dict = {}
   with open( ref_file, 'r') as mr:      
      cur_chr = ''
      cur_seq = [];
      while True:
         line = mr.readline();
         if not line: break;
         line = line.strip();
         if len(line)==0: continue;

         if line[0]=='>':
             if len(cur_seq)>0:
                 if '-' in cur_chr or '_' in cur_chr:
                    if len(sys.argv)>3 and sys.argv[3] in ['1', 'y', 'yes', 'Y', 'YES']:
                       ref_dict[ cur_chr ] = ''.join( cur_seq );
                 else:
                    ref_dict[ cur_chr ] = ''.join( cur_seq );
                 #if cur_chr not in orderedchr: orderedchr.append( cur_chr )
                 #if cur_chr not in orderedchrnum: orderedchrnum.append( cur_chr )
             cur_chr = line.split()[0][1:]
             cur_seq = []
         else:
             cur_seq.append(line.upper());

      if len(cur_seq)>0:
         if '-' in cur_chr or '_' in cur_chr:
            if len(sys.argv)>3 and sys.argv[3] in ['1', 'y', 'yes', 'Y', 'YES']:
               ref_dict[ cur_chr ] = ''.join( cur_seq );
         else:
            ref_dict[ cur_chr ] = ''.join( cur_seq );
         #if cur_chr not in orderedchr: orderedchr.append( cur_chr )
         #if cur_chr not in orderedchrnum: orderedchrnum.append( cur_chr )

   in_list = [0, 0]
   for _rk in ref_dict:
      if _rk in orderedchr: in_list[0] = in_list[0]+1
      if _rk in orderedchrnum: in_list[1] = in_list[1]+1
   
   for_sav_chrs = orderedchr if in_list[0]>in_list[1] else orderedchrnum
   print( len(ref_dict), in_list );
   for _rk in ref_dict:
      if _rk not in for_sav_chrs: for_sav_chrs.append( _rk );
   cpg_id = 0;
   mw = open( op_file, 'w')
   for _rk in for_sav_chrs:
      if _rk in ref_dict:
         cpg_id = write_cpg( mw, _rk, ref_dict[_rk], cpg_id);       
   mw.close();

