class KSimplex:
    def __init__(self, listofvertices):
        self.kvertices = listofvertices
        self.k = len(listofvertices) - 1
        self.name = str(self.k) + '-simplex: ' + str(listofvertices)
        self.id = -1  # This id is used as index while building transformation matrix

    def __str__(self):
        # return self.name
        return str(self.kvertices)

    def __eq__(self, other):
        '''
        Check whether two k-simplex are the same irrespective of the orientation
        :param other: Another KSimplex
        :return: true if both simplex are same or may be just different orientation
        '''
        if self.k == other.k:
            return sorted(self.kvertices) == sorted(other.kvertices)
        return False


class SimplicialComplex:
    def __init__(self):
        self.simplex = []  # Stores all the simplices in the complex
        self.tableofksimplex = {}  # key = k , value = list of k-simplices in the simplicial_complex
        self.maxK = 0  # Keep track of highest Dimensional simplex in the complex
        self.count_id = {}  # for assigning unique id to each set of simplex

    def get_allkth_simplices(self, k):
        return self.tableofksimplex.get(k, [])

    def add_simplex_fromfile(self, filename):
        '''
        Takes filename as input stream for simplices
        :param filename:
        :return:
        '''
        with open(filename, 'r') as fp:
            while 1:
                line = fp.readline()
                if not line:
                    break
                ksimplex_obj = KSimplex([int(v) for v in line.split()])  # Building the K-simplex object
                self.add_simplex(ksimplex_obj)

    def add_simplex(self, ksimplex):
        '''
        Add a k-simplex to the simplicial complex. Can be added in any order and any dimension.
        :param ksimplex:
        :return: None
        '''
        assert isinstance(ksimplex, KSimplex)
        self.simplex.append(ksimplex)
        if self.tableofksimplex.get(ksimplex.k, None) is None:
            self.tableofksimplex[ksimplex.k] = []
            self.count_id[ksimplex.k] = 0
        self.tableofksimplex[ksimplex.k].append(ksimplex)
        ksimplex.id = self.count_id[ksimplex.k]  # assigning id
        self.count_id[ksimplex.k] += 1  # increasing id for the next k-simplex's id to be id+1

        self.maxK = [self.maxK, ksimplex.k][
            ksimplex.k > self.maxK]  # update maxK if a higher dimensional simplex is added

    def __str__(self):
        toreturn = ''
        for k, simplices in self.tableofksimplex.items():
            toreturn += str(k) + ': ' + str([str(x) for x in simplices]) + ' '
            toreturn += '\n'
        return toreturn
