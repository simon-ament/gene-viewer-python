from abc import ABC, abstractmethod


class Processor(ABC):
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
        for transcript_id, transcript_regions in data["regions"].items():
            exon_junctions = list(
                filter(lambda x: x["type"] == "exon_junction", transcript_regions)
            )

            if not exon_junctions:
                continue

            sorted_exon_junctions = sorted(exon_junctions, key=lambda x: x["start"])

            merged_exon_junctions = []
            last_exon_junction = sorted_exon_junctions[0]
            for exon_junction in sorted_exon_junctions[1:]:
                if exon_junction["exon_number"] == last_exon_junction["exon_number"]:
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
                list(filter(lambda x: x["type"] != "exon_junction", transcript_regions))
                + merged_exon_junctions
            )

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
        exons = []
        for transcript_data in data["regions"].values():
            exons.extend(list(filter(lambda x: x["type"] == "exon", transcript_data)))

        sorted_exons = sorted(exons, key=lambda x: x["start"])
        sorted_sequences = sorted(
            data["sequences"], key=lambda x: (x["start"], len(x["sequence"]))
        )

        # iterate through exons and sequences in parallel to restrict sequences to only exons
        exon_idx = 0
        sequence_idx = 0
        while exon_idx < len(sorted_exons) and sequence_idx < len(sorted_sequences):
            exon = sorted_exons[exon_idx]
            sequence = sorted_sequences[sequence_idx]

            exon_start = exon["start"]
            exon_end = exon["end"]
            seq_start = sequence["start"]
            seq_end = seq_start + len(sequence["sequence"]) - 1

            if seq_end < exon_start:
                # sequence is before the exon, move to the next sequence
                sequence_idx += 1
            elif seq_start > exon_end:
                # sequence is after the exon, move to the next exon
                exon_idx += 1
            else:
                # sequence overlaps with the exon, restrict it to the exon boundaries
                restricted_start = max(seq_start, exon_start)
                restricted_end = min(seq_end, exon_end)
                restricted_sequence = sequence["sequence"][
                    restricted_start - seq_start : restricted_end - seq_start + 1
                ]
                data["sequences"].append(
                    {"start": restricted_start, "sequence": restricted_sequence}
                )
                if restricted_end == seq_end:
                    # sequence is fully contained within the exon, move to the next sequence
                    sequence_idx += 1
                else:
                    # sequence extends beyond the exon, move to the next exon
                    exon_idx += 1

        return data


# Other ideas
# - UTR and CDS to exons
