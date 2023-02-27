Created 2/22/23
author: etmu9498

Purpose of these scripts:
- to process and save 2022 crl data in a much simpler format, where it can be
  ingested into statistical calculations
- previous scripts have lots of redundant code that relies on the presence of
  tdr and flight level datasets for processing. These axes will be optional when
  processing new datasets!


save_crl_data.py script goals:
- try not to rely on manual metadata inputs
- save datasets as "[date][p-3]_[name]_processed.nc"
