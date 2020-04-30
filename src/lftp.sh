# commands for downloading from FTPS with lftp

cd /secondary/projects/bbc/research/PFEG_20200427_CDseq
screen -S lftp
lftp ftp://pinkpig96:65VtdA2L@ftps.fulgentgenetics.com:21
set ssl:verify-certificate false
mirror FT-SE6114
