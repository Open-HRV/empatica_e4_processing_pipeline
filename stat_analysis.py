import bioread
import os.path
import shutil
from os import listdir

import pandas as pd

dir = os.getcwd() + r"\results\\"

results = [name for name in listdir(dir)]

print(results)
agg = pd.DataFrame();
for result in results:
    hrv = pd.read_csv(dir + r"\\" + result, sep=",")
    print(result)
    agg[result[4:-4] + "_HRV_RMSSD"] = hrv["HRV_RMSSD"]
    agg[result[4:-4] + "_HRV_MeanNN"] = hrv["HRV_MeanNN"]
    agg[result[4:-4] + "_HRV_SDNN"] = hrv["HRV_SDNN"]
    agg[result[4:-4] + "_HRV_MedianNN"] = hrv["HRV_MedianNN"]
    agg[result[4:-4] + "_HRV_LF"] = hrv["HRV_LF"]
    agg[result[4:-4] + "_HRV_HF"] = hrv["HRV_HF"]
    agg[result[4:-4] + "_HRV_VHF"] = hrv["HRV_VHF"]
    agg[result[4:-4] + "_HRV_LFHF"] = hrv["HRV_LFHF"]
    agg[result[4:-4] + "_HRV_ApEn"] = hrv["HRV_ApEn"]
    agg[result[4:-4] + "_HRV_SampEn"] = hrv["HRV_SampEn"]

agg.to_csv("biopac_vs_e4_full.csv", mode='w', header=True, index=True)



