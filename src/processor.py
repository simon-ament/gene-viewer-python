from abc import ABC, abstractmethod

class Processor(ABC):
    def __init__(self, ):
        pass
        
    @property
    @abstractmethod
    def id(self):
        pass

    @property
    @abstractmethod
    def version(self):
        pass

    @property
    @abstractmethod
    def apply_to(self):
        # which data types this processor will be applied to ("regions", "sequences", "tracks", or "probes")
        # will be applied in the specified order
        pass

    @property
    @abstractmethod
    def dependencies(self):
        # a dictionary specifying for each data type, which other data types this processor depends on
        # e.g. {"regions": ["sequences"], "tracks": ["regions", "sequences"]}
        pass

    @abstractmethod
    def process(self, data):
        pass

class ProcessorExonJunctions(Processor):
    id = "exon_junctions"
    version = "1"
    apply_to = ["regions"]
    dependencies = {}

    def __init__(self, ):
        super().__init__()

    def process(self, data):
        # merge exon junctions with same exon_number into single exon
        return data

class ProcessorIntronGaps(Processor):
    id = "intron_gaps"
    version = "1"
    apply_to = ["regions"]
    dependencies = {}

    def __init__(self, ):
        super().__init__()

    def process(self, data):
        # add introns in gaps between exons
        return data

class ProcessorExonSequencesOnly(Processor):
    id = "exon_sequences_only"
    version = "1"
    apply_to = ["regions"]
    dependencies = {"regions": ["sequences"]}

    def __init__(self, ):
        super().__init__()

    def process(self, data):
        # restrict sequences to only exons
        return data
