from abc import ABC, abstractmethod


class Processor(ABC):
    @property
    @abstractmethod
    def id(self):
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
    id = "exon_junctions_v1"
    input = ("regions",)
    output = ("regions",)

    def __init__(
        self,
    ):
        pass

    def process(self, data):
        # merge exon junctions with same exon_number into single exon
        for transcript_id, transcript_regions in data["regions"].items():
            exon_junctions = list(
                filter(lambda x: x["type"] == "exonexonjunction", transcript_regions)
            )

            if not exon_junctions:
                continue

            sorted_exon_junctions = sorted(exon_junctions, key=lambda x: x["start"])

            merged_exon_junctions = []
            last_exon_junction = sorted_exon_junctions[0]
            for exon_junction in sorted_exon_junctions[1:]:
                exon_number = exon_junction["exon_number"] # first exon number
                last_exon_number = last_exon_junction["exon_number"] # second exon number
                if exon_number == last_exon_number:
                    # merge exon junctions with same exon_number
                    last_exon_junction["end"] = max(
                        last_exon_junction["end"], exon_junction["end"]
                    )
                    last_exon_junction["type"] = "exon"
                else:
                    merged_exon_junctions.append(last_exon_junction)
                    last_exon_junction = exon_junction
            merged_exon_junctions.append(last_exon_junction)

            # replace exon junctions with merged exon junctions
            data["regions"][transcript_id] = (
                list(filter(lambda x: x["type"] != "exonexonjunction", transcript_regions))
                + merged_exon_junctions
            )

        return data


class ProcessorIntronGaps(Processor):
    id = "intron_gaps_v1"
    input = ("regions",)
    output = ("regions",)

    def __init__(
        self,
    ):
        pass

    def process(self, data):
        # add introns in gaps between exons
        for transcript_id, transcript_regions in data["regions"].items():
            exons = list(filter(lambda x: x["type"] == "exon", transcript_regions))
            if not exons:
                continue

            # Sort exons by start position
            sorted_exons = sorted(exons, key=lambda x: x["start"])

            introns = []
            for i in range(len(sorted_exons) - 1):
                intron_start = sorted_exons[i]["end"] + 1
                intron_end = sorted_exons[i + 1]["start"] - 1
                if intron_start <= intron_end:
                    introns.append(
                        {"start": intron_start, "end": intron_end, "type": "intron"}
                    )

            # Add introns to the transcript regions
            data["regions"][transcript_id].extend(introns)

        return data


class ProcessorExonSequencesOnly(Processor):
    id = "exon_sequences_only_v1"
    input = ("regions", "sequences")
    output = ("sequences",)

    def __init__(
        self,
    ):
        pass

    def process(self, data):
        # restrict sequences to only exons
        exons = []
        for transcript_data in data["regions"].values():
            exons.extend(list(filter(lambda x: x["type"] == "exon", transcript_data)))

        sorted_exons = sorted(exons, key=lambda x: x["start"])
        sorted_sequences = sorted(
            data["sequences"], key=lambda x: (x["start"], len(x["sequence"]))
        )

        # iterate through sequences and exons in parallel to restrict sequences to only exons
        exon_idx = 0
        sequence_idx = 0
        selection_start = 0
        selection_end = 0
        sequences = []

        while sequence_idx < len(sorted_sequences) and exon_idx < len(sorted_exons):
            # 1. find next selection start, i.e. the next position where both a sequence and an exon exist
            while (
                sequence_idx < len(sorted_sequences)
                and exon_idx < len(sorted_exons)
                and (
                    sorted_sequences[sequence_idx]["start"]
                    > sorted_exons[exon_idx]["end"]
                    or sorted_exons[exon_idx]["start"]
                    > sorted_sequences[sequence_idx]["start"] + len(sorted_sequences[sequence_idx]["sequence"]) - 1
                )
            ):
                if sorted_sequences[sequence_idx]["start"] < sorted_exons[exon_idx]["start"]:
                    sequence_idx += 1
                else:
                    exon_idx += 1

            selection_start = max(
                sorted_sequences[sequence_idx]["start"], sorted_exons[exon_idx]["start"]
            )

            # 2. find next selection end, i.e. the next position where the sequence ends or no continuous exon exists
            exon_end = sorted_exons[exon_idx]["end"]
            while (
                exon_idx + 1 < len(sorted_exons)
                and sorted_exons[exon_idx + 1]["start"] <= exon_end + 1
            ):
                exon_idx += 1
                exon_end = sorted_exons[exon_idx]["end"]

            selection_end = min(
                sorted_sequences[sequence_idx]["start"] + len(sorted_sequences[sequence_idx]["sequence"]) - 1, exon_end
            )
            
            # 3. select the sequence from selection_start to selection_end and add it to the sequences list
            if selection_start <= selection_end:
                sequence = sorted_sequences[sequence_idx]["sequence"][
                    selection_start - sorted_sequences[sequence_idx]["start"] : selection_end
                    - sorted_sequences[sequence_idx]["start"]
                    + 1
                ]
                sequences.append({"start": selection_start, "sequence": sequence})

            # 4. move to the next sequence or exon
            if sorted_exons[exon_idx]["end"] <= selection_end:
                exon_idx += 1
            else:
                sequence_idx += 1

        data["sequences"] = sequences
        return data


# Other ideas
# - UTR and CDS to exons
