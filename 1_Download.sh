for i in GSM5678319_BM-UPN01 GSM5678320_BM-UPN02 GSM5678321_BM-UPN03 GSM5678326_CD34-Old4 GSM5678327_CD34-Old5; do
  echo "Processing: $i"

  gsm_id="${i%%_*}"
  echo "$gsm_id"
  echo "$i"
  wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5678nnn/${gsm_id}/suppl/${i}_counts.tar.gz
