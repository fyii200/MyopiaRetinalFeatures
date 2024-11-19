import os
import pandas as pd
from OD.dataset import DatasetGenerator
from OD.model import trainedModel
from torch.utils.data import DataLoader
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore")
join = os.path.join

# Read dataframe
df = pd.read_csv("cleaned_data_long_PM_cohort.csv")

# Initialise UK Biobank dataset & dataloader
UKBdataset = DatasetGenerator(data_frame = df)
UKBdataloader = test_loader = DataLoader(dataset=UKBdataset, batch_size=1, shuffle=False, num_workers=8, pin_memory=False)

# Load the trained model (DeepLabV3 with a MobileNet-Large backbone)
model = trainedModel(savedMaskDir = "segmentedOD")
model.load(dir_path = "OD")

# Start segmentation (predicted binary masks are saved automatically)
for (image, image_name) in tqdm(UKBdataloader):
    # Note: image_name is a tuple containing an image name,
    # so need to index it with [0] to extract the string.
    predictedMask = model.segment(image, image_name[0])