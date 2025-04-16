import pandas as pd
import numpy as np

import pandas as pd
import numpy as np

def file_reader(path, column=0, base=16):

    print(f"Reading: {path}")
    try:
        df = pd.read_csv(path, header=None, usecols=[0, 1])
        df.columns = ["h", "l"]
        data = np.vectorize(int)(df["h"].to_numpy(), base) , np.vectorize(int)(df["l"].to_numpy(), base)
    except:
        print("Error with the reading")
    return data
