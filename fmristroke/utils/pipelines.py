import json

class Pipeline:
    def __init__(self, pipeline = None, confounds = None, demean = None, clean_spec = None, desc=None):
        self.pipeline = pipeline
        self.confounds_spec = confounds
        self.demean = demean
        self.clean_specs = clean_spec
        self._desc = desc
        
    @classmethod
    def from_json(cls, json_file):
        with open(json_file, "rb") as f:
            pipe = json.load(f)
        return(cls(**pipe))
    def __str__(self):
        return self.pipeline
    

class Pipelines:
    @staticmethod
    def check_pipeline(pipeline):
        if isinstance(pipeline, Pipeline):
            return pipeline
        
        try:
            pipeline = Pipeline.from_json(pipeline)
        except Error:
            return Pipeline()
        
    def __init__(self, pipelines=None):
        self._pipelines = pipelines
    @property
    def pipelines(self):
        return self._pipelines
    def add(self, val):
        if val not in self._pipelines:
            self._pipelines += [self.check_pipeline(val)]
    def get_pipelines(self):
        out = []
        for s in self.pipelines:
            out.append(s.pipeline)
        return out
    def get_confounds_spec(self):
        out = []
        for s in self.pipelines:
            out.append(s.confounds_spec)
        return out
    def get_demean(self):
        out = []
        for s in self.pipelines:
            out.append(s.demean)
        return out
    def get_clean_specs(self):
        out = []
        for s in self.pipelines:
            out.append(s.clean_specs)
        return out
    def __str__(self):
        pipelines = ", ".join([str(s) for s in self.pipelines]) or "<none>."
        return "Pipelines: %s" % pipelines
    