import os,sys
import torch

from models_instanceNorm_reco_2chan_quadTask import ResBlock, ResNet34

def larpid_loadmodel(model_path, device):
    print("LOADING RESNET MODEL")
    model = ResNet34(2, ResBlock, outputs=5)
    if device == "cpu":
        checkpoint = torch.load(model_path, map_location=torch.device('cpu'))
    else:
        checkpoint = torch.load(model_path)
        
    try:
        model.load_state_dict(checkpoint['model_state_dict'])
    except:
        model.module.load_state_dict(checkpoint['model_state_dict'])
    model.to(args.device)
    model.eval()    
    return model
