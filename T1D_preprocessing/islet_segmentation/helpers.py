import cv2
import numpy as np
import pandas as pd
import torch

from pathlib import Path
from tqdm import tqdm

import albumentations as A
from torch import nn
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
from torch.optim import lr_scheduler
from torchmetrics.classification import BinaryF1Score, BinaryAccuracy
import segmentation_models_pytorch as smp
from segmentation_models_pytorch.losses import DiceLoss



# Data augmentation functions
def get_train_augmentations(img_size):
    return A.Compose([
        A.Resize(img_size, img_size),
        A.VerticalFlip(p = 0.5),
        A.HorizontalFlip(p = 0.5),
        A.RandomRotate90(p = 0.5),
        # A.OneOf([
        #     A.ElasticTransform(alpha=120, sigma=120 * 0.05, alpha_affine=120 * 0.03, p=0.2),
        #     A.GridDistortion(p=0.2),
        #     A.OpticalDistortion(distort_limit=2, shift_limit=0.2, p=1)                  
        #     ], p=0.5),
        A.RandomBrightnessContrast(p = 0.5)
    ])

def get_test_augmentations(img_size):
    return A.Compose([
        A.Resize(img_size, img_size)
    ])


# Islet dataset class
class IsletDataset(Dataset):
    
    def __init__(self, df, augmentations):
        self.df = df
        self.augmentations = augmentations
        
    def __len__(self):
        return len(self.df)
    
    def __getitem__(self, idx):
        
        row = self.df.iloc[idx]
        image_path = row.dir_images / row.images
        mask_path = row.dir_masks / row.masks
        
        image = cv2.imread(str(image_path), cv2.IMREAD_UNCHANGED)
        image = np.expand_dims(image, axis = -1)
        
        mask = cv2.imread(str(mask_path), cv2.IMREAD_UNCHANGED)
        mask = np.expand_dims(mask, axis = -1)
        
        if self.augmentations:
            data = self.augmentations(image = image, mask = mask)
            image = data['image']
            mask = data['mask']

        image = np.transpose(image, (2,0,1)).astype(np.float32)
        mask = np.transpose(mask, (2,0,1)).astype(np.float32)
        
        image = torch.Tensor(image) / 255.0
        mask = torch.round(torch.Tensor(mask))
        
        return image, mask
    
class PredictDataset(Dataset):
    
    def __init__(self, df, augmentations):
        self.df = df
        self.augmentations = augmentations
        
    def __len__(self):
        return len(self.df)
    
    def __getitem__(self, idx):
        
        row = self.df.iloc[idx]
        image_path = row.dir_images / row.images        
        image = cv2.imread(str(image_path), cv2.IMREAD_UNCHANGED)
        image = np.expand_dims(image, axis = -1)
        
        if self.augmentations:
            data = self.augmentations(image = image)
            image = data['image']

        image = np.transpose(image, (2,0,1)).astype(np.float32)
        image = torch.Tensor(image) / 255.0

        return image


# Training and evaluation functions
def training_function(data_loader, model, optimizer, device):
    
    model.train()
    total_loss = 0.0
    
    for images, masks in tqdm(data_loader):
        images = images.to(device)
        masks = masks.to(device)
        
        optimizer.zero_grad()
        logits, loss = model(images, masks)
        loss.backward()
        optimizer.step()
        
        total_loss += loss.item()
        
    return total_loss / len(data_loader)

def evaluation_function(data_loader, model, device):
    
    model.eval()
    total_loss = 0.0
    total_f1 = 0.0
    total_acc = 0.0
    # Push metrics to device.
    binary_acc = BinaryAccuracy().to(device)
    binary_f1 = BinaryF1Score().to(device)

    with torch.no_grad():
        for images, masks in tqdm(data_loader):
            images = images.to(device)
            masks = masks.to(device)

            # Loss per batch
            logits, loss = model(images, masks)
            total_loss += loss.item()
            
            # Accuracy per batch.
            mask_pred = torch.sigmoid(logits)
            mask_pred = (mask_pred > 0.5) * 1.0
            total_acc += binary_acc(preds = mask_pred, target = masks)

            # F1 score per batch.
            total_f1 += binary_f1(preds = mask_pred, target = masks)
    
    # return loss and accuracy per Epoch.
    return total_loss / len(data_loader), total_acc / len(data_loader), total_f1 / len(data_loader)


# Islet segmentation model
class IsletModel(nn.Module):
    
    def __init__(self):
        super(IsletModel, self).__init__()
        
        # Architecture (arc)
        self.arc = smp.Unet(
            encoder_name = "timm-efficientnet-b0",
            encoder_weights = "imagenet",
            in_channels = 1,  # number of channel in input images
            classes = 1,
            activation = None
        )
        
    def forward(self, images, masks = None):
        
        logits = self.arc(images)
        
        if masks != None:
            loss1 = DiceLoss(mode = 'binary')(logits, masks)
            loss2 = nn.BCEWithLogitsLoss()(logits, masks)
            return logits, loss1 + loss2
        
        return logits