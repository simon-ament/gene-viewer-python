from abc import ABC, abstractmethod

class IndexedFile(ABC):
    def __init__(self, file_path):
        self.file_path = file_path
        self.file = open(file_path, 'r')
        self.index = self.build_index()

    @abstractmethod
    def build_index(self):
        pass

    @abstractmethod
    def get(self, key):
        pass

    @abstractmethod
    def keys(self):
        pass

    def __del__(self):
        self.file.close()

# key: SeqID, return: (header, sequence)
class IndexedFastaFile(IndexedFile):
    def build_index(self):
        index = {}
        position = 0
        for line in self.file:
            if line.startswith('>'):
                key = line[1:].strip().split()[0]  # Use the first word after '>' as the key
                index[key] = position
            position += len(line)
        return index

    def get(self, key):
        if key not in self.index:
            raise KeyError(f"Key '{key}' not found in index.")
        position = self.index[key]
        self.file.seek(position)
        header = self.file.readline().strip()
        sequence_lines = []
        while True:
            line = self.file.readline()
            if not line or line.startswith('>'):
                break
            sequence_lines.append(line.strip())
        return header, ''.join(sequence_lines)

    def keys(self):
        return list(self.index.keys())

# key: gene ID, return: list of features
class IndexedGTFFile(IndexedFile):
    def build_index(self, genes):
        pass

    def get(self, key):
        pass

    def keys(self):
        return list(self.index.keys())

# key: gene ID, return: list of features
class IndexedBEDFile(IndexedFile):
    def build_index(self, genes):
        pass

    def get(self, key):
        pass

    def keys(self):
        return list(self.index.keys())
