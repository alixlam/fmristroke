import json
import os

class Atlas:
    def __init__(self, atlas = None, space = None, labels = None, mask_file=None, desc=None):
        self.atlas = atlas
        self.space = space
        self.labels = labels
        self.mask_file = mask_file
        self._desc = desc
        
    @classmethod
    def from_json(cls, json_file):
        with open(json_file, "rb") as f:
            pipe = json.load(f)
        # Mask file path 
        mask_file_path = os.path.join(os.path.dirname(json_file), pipe["mask_file"])
        pipe["mask_file"] = mask_file_path
        return(cls(**pipe))
    def __str__(self):
        return self.atlas
    
class Atlases:
    @staticmethod
    def check_atlas(atlas):
        if isinstance(atlas, Atlases):
            return atlas
        
        try:
            atlas = Atlas.from_json(atlas)
        except Error:
            return Atlas()
        
    def __init__(self, atlases=None):
        self._atlases = atlases
    @property
    def atlases(self):
        return self._atlases
    def add(self, val):
        if val not in self._atlases:
            self._atlases += [self.check_atlas(val)]
    def get_atlases(self):
        out = []
        for s in self.atlases:
            out.append(s.atlas)
        return out
    def get_labels(self):
        out = []
        for s in self.atlases:
            out.append(s.labels)
        return out
    def get_spaces(self):
        out = []
        for s in self.space:
            out.append(s.demean)
        return out
    def get_mask_files(self):
        out = []
        for s in self.atlases:
            out.append(s.mask_file)
        return out
    def __str__(self):
        atlases = ", ".join([str(s) for s in self.atlases]) or "<none>."
        return "Atlases: %s" % atlases
    