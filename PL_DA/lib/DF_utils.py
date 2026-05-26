
import os, sys, yaml


def _yaml2dict(cfgFile) :
    with open(cfgFile) as fstream :
        cfgDict = yaml.safe_load(fstream)
    
    return cfgDict

