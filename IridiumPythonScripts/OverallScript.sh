#!/bin/bash

set -e

python3 IridiumGeomagneticFieldDataMaps.py

python3 AccioSV.py

python3 SmoothSV_AccioJerk.py

echo 'I am done smoothing'

cd ./JerkData/

# rm ./MasterJerk*.mat # Uncomment if you are running the script not for the first time to remove existing files

sed -i 's/JerkTimes_Bphi/JerkTimes_Br/' JoinJerkData.py
sed -i 's/MasterJerkValley_Bphi.mat/MasterJerkValley_Br.mat/' JoinJerkData.py
sed -i 's/MasterJerkPeak_Bphi.mat/MasterJerkPeak_Br.mat/' JoinJerkData.py

python3 JoinJerkData.py 

echo 'R component joined'

sed -i 's/JerkTimes_Br/JerkTimes_Btheta/' JoinJerkData.py 
sed -i 's/MasterJerkValley_Br.mat/MasterJerkValley_Btheta.mat/' JoinJerkData.py
sed -i 's/MasterJerkPeak_Br.mat/MasterJerkPeak_Btheta.mat/' JoinJerkData.py

python3 JoinJerkData.py 

echo 'Theta component joined'

sed -i 's/JerkTimes_Btheta/JerkTimes_Bphi/' JoinJerkData.py 
sed -i 's/MasterJerkValley_Btheta.mat/MasterJerkValley_Bphi.mat/' JoinJerkData.py
sed -i 's/MasterJerkPeak_Btheta.mat/MasterJerkPeak_Bphi.mat/' JoinJerkData.py

python3 JoinJerkData.py

echo 'Phi component joined'

sed -i 's/_Btheta/_Br/' JoinPeakValley.py

python3 JoinPeakValley.py

echo 'Joined Peak & Valley for R'

sed -i 's/_Br/_Bphi/' JoinPeakValley.py

python3 JoinPeakValley.py

echo 'Joined Peak & Valley for Phi'

sed -i 's/_Bphi/_Btheta/' JoinPeakValley.py

python3 JoinPeakValley.py

echo 'Joined Peak & Valley for Theta'