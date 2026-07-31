from abc import ABC, abstractmethod


class Processor(ABC):
    def __init__(
        self,
    ):
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
    def input(self):
        pass

    @property
    @abstractmethod
    def output(self):
        pass

    @abstractmethod
    def process(self, data):
        pass


class ProcessorExonJunctions(Processor):
    id = "exon_junctions"
    version = "1"
    input = ("regions",)
    output = ("regions",)

    def __init__(
        self,
    ):
        pass

    def process(self, data):
        # merge exon junctions with same exon_number into single exon
        return data


class ProcessorIntronGaps(Processor):
    id = "intron_gaps"
    version = "1"
    input = ("regions",)
    output = ("regions",)

    def __init__(
        self,
    ):
        pass

    def process(self, data):
        # add introns in gaps between exons
        return data


class ProcessorExonSequencesOnly(Processor):
    id = "exon_sequences_only"
    version = "1"
    input = ("regions", "sequences")
    output = ("regions",)

    def __init__(
        self,
    ):
        pass

    def process(self, data):
        # restrict sequences to only exons
        return data
