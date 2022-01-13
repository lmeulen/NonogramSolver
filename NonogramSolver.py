import time
import cProfile


class Nonogram:
    """
    This class implements a solver for Nonogram puzzels.
    Typical usage:
        n = Nonogram()
        n.initialize_from_file('puzzle.txt', verbose=True)
        n.solve()
        n.print(fields=False, print_board=True)
    """
    no_rows = None
    no_cols = None
    board = None
    col_hints = None
    row_hints = None
    possible_rows = None
    flipped = False

    def initialize(self, row_count, col_count):
        """
        Initialize a Nonogram. Only square nonograms are supported.
        :param row_count: Number of  rows in the nonogram
        :param col_count: Number of columns in the nonogram
        :return:
        """
        self.no_cols = col_count
        self.no_rows = row_count
        self.col_hints = None
        self.row_hints = None
        self.possible_rows = None
        self.flipped = False
        self.board = []
        for rowid in range(0, self.no_rows):
            self.board.append([0] * self.no_cols)

    def initialize_from_file(self, filename, verbose=False):
        """
        Initialize a nonogram from a file.
        File format:
            - first row is the number of rows
            - second row is the number of columns
            - then for each row the hints, seperated by space
            - then for each column the hints, seperated by spaces.
            Total number of lines is (2 + no_rows + no_columns)
        Prints the initialised nonogram.
        :param filename: Filename of file with nonogram puzzle specification
        :param verbose: If True, print the initialised class
        :return:
        """
        with open(filename) as f:
            # Get sizes
            row_count = int(f.readline())
            col_count = int(f.readline())
            self.initialize(row_count, col_count)
            # Read row hints, one row per line
            self.row_hints = []
            for i in range(0, self.no_rows):
                line = f.readline()
                numbers = [int(i) for i in line.split()]
                self.row_hints.append(numbers)
            # Read column hints, one column per line
            self.col_hints = []
            for i in range(0, self.no_cols):
                line = f.readline()
                numbers = [int(i) for i in line.split()]
                self.col_hints.append(numbers)
        if verbose:
            self.print(fields=True, print_board=False)

    def print(self, fields=False, print_board=True):
        """
        Print the contents of the class
        :param fields: Print all relevant class fields (except board, controlled by other parameter)
        :param print_board: Print the board stored in the class instance
        :return:
        """
        if fields:
            print('Size          : ' + str(self.no_rows) + ' X ' + str(self.no_cols))
            print('Column hints  : ' + str(self.col_hints))
            print('Row hints     : ' + str(self.row_hints))
        if fields and print_board:
            print('Board')
        if print_board and self.board:
            chars = ['.', 'X']
            for r in self.board:
                for c in r:
                    print(chars[c], end=' ')
                print()

    def set_row(self, rownumber, values):
        """
        Set the board values for the specfied rows.
        Rows start counting at 0.
        :param rownumber: Rownumber to set
        :param values: Values for the row. A list of integers with length = size
        :return:
        """
        self.board[rownumber] = values.copy()

    def get_board(self):
        """
        Return a copy of the stored board as array of array of ints
        Hereby 0 represents an empty cell and 1 represents a filled cell
        :return: Array of array of ints
        """
        return self.board.copy()

    def get_size(self):
        """
        Return the size of the board
        The size equals the number of rows and the number of columns
        :return: The size of the board as integer
        """
        return self.no_rows, self.no_cols

    def set_col_hints(self, hints):
        """
        Set the column hints. The hinsts are given in a list of columns
        where each column is a list of numbers. E.g. [[1,2],[2]]
        :param hints: Hints as list of list of ints
        :return:
        """
        self.col_hints = hints.copy()

    def set_row_hints(self, hints):
        """
        Set the row hints. The hinsts are given in a list of rows
        where each row is a list of numbers. E.g. [[1,2],[2]]
        :param hints: Hints as list of list of ints
        :return:
        """
        self.row_hints = hints.copy()

    def get_row_descriptor(self, rowid):
        """
        Return the row descriptor (hints) of a given row of the
        board stored in the class. Return as a list of numbers.
        E.g. ..XX..X.. reswults in [2, 1]
        :param rowid: Row ID, row ID's start at 0
        :return: list of integers describigng the row
        """
        descriptor = []
        previous_field = 0
        count = 0
        '''
        We iterate over the row. 
        If a '1' is found
        - if previous element is a zero, a new serie of '1' is found,
        - if previous element is a one, the current serie is extende with 1
        If a '0' is found, the current serie is added to the descrptor
        If at the end of the row a serie is in progress, it is added to the list
        '''
        for colid in range(0, self.no_cols):
            if self.board[rowid][colid] == 1:
                if previous_field < 1:
                    count = 1
                else:
                    count += 1
            elif self.board[rowid][colid] == 0 and previous_field == 1:
                descriptor.append(count)
                count = 0
            previous_field = self.board[rowid][colid]
        if count > 0:
            descriptor.append(count)
        return descriptor

    def get_col_descriptor(self, colid, no_rows=None):
        """
        Return the column descriptor (hints) of a given column of the
        board stored in the class. Return as a list of numbers.
        The number of rows used can be limited.
        E.g. ..XX..X.. (transposed)reswults in [2, 1]
        :param colid: Column ID, ID's start at 0
        :param no_rows: Number of rows to use
        :return: list of integers describigng the column
        """
        descriptor = []
        previous_field = 0
        count = 0
        if not no_rows:
            no_rows = self.no_rows
        '''
        We iterate over the column. 
        If a '1' is found
        - if previous element is a zero, a new serie of '1' is found,
        - if previous element is a one, the current serie is extende with 1
        If a '0' is found, the current serie is added to the descrptor
        If at the end of the row a serie is in progress, it is added to the list
        '''
        for rowid in range(0, no_rows):
            current_field = self.board[rowid][colid]
            if current_field == 1:
                if previous_field < 1:
                    count = 1
                else:
                    count += 1
            elif previous_field == 1:  # and current_field = 0
                descriptor.append(count)
                count = 0
            previous_field = current_field
        if count > 0:
            descriptor.append(count)
        return descriptor

    def get_row_descriptors(self):
        """
        Returns a list with the row descriptor for all rows
        E.g.
        . . X
        X . X
        . X X
        returns [[1], [1,1], [2]]
        :return: List of lists of numbers, list of row descriptors
        """
        res = []
        for rowid in range(0, self.no_rows):
            res.append(self.get_row_descriptor(rowid))
        return res

    def get_col_descriptors(self):
        """
        Returns a list with the column descriptor for all columns
        E.g.
        . . X
        X . X
        . X X
        returns [[1], [1], [3]]
        :return: List of lists of numbers, list of column descriptors
        """
        res = []
        for colid in range(0, self.no_cols):
            res.append(self.get_col_descriptor(colid))
        return res

    def __gen_row_next_part(self, row, hints_index, hints, spaces, xtra_spaces_limit):
        """
        Recursively create possible rows
        :param row: Row so far
        :param hints_index: Current position in the description
        :param hints: Row description, filled fields
        :param spaces: Spaces description
        :param xtra_spaces_limit: Additional spaces possible
        :return:
        """
        result = []
        if hints_index < len(hints):
            # Add all possible number of spaces at the current position
            for additional_spaces in range(spaces[hints_index], spaces[hints_index] + xtra_spaces_limit + 1):
                # Extended row with spaces and next hint
                next_row = row + [0] * additional_spaces + [1] * hints[hints_index]
                # Recursively add next elements
                sub_res = self.__gen_row_next_part(next_row, hints_index + 1, hints, spaces, xtra_spaces_limit)
                for r in sub_res:
                    result.append(r)
        elif len(row) <= self.no_cols:
            # All hints already added, time to return result
            # Return possible solution, after filling it with 0's to match size
            return [(row + [0] * (self.no_cols - len(row)))]
        return result

    def get_row_possibilities(self, rowid):
        """
        Generate all possible rows for a given row.
        :param rowid: Row ID (starting with 0
        :return: List of all possible rows conform the row hints
        """
        hints = self.row_hints[rowid]
        # Create a list that stores the number of empty fields before each serie of filled fields
        # We start with one space between each serie of 1's, except for the first serie where we start
        # with zero empty fields.
        spaces = [1] * len(hints)
        spaces[0] = 0
        # Calculate the number of possible empty fields to add.
        xtra_spaces = self.no_cols - sum(hints) - len(hints) + 1
        possibilities = self.__gen_row_next_part([], 0, hints, spaces, xtra_spaces)
        return possibilities

    def generate_possible_rows(self):
        """
        Generate all possible rows for the given row hints
        Stored as a list with an item per row. This item contains a list of
        all possible rows. Stored in the class and returned.
        :return: The list of possible rows per row
        """
        self.possible_rows = []
        for rowid in range(0, self.no_rows):
            self.possible_rows.append(self.get_row_possibilities(rowid))
        return self.possible_rows

    def check_board_vs_col_constraints(self):
        """
        Check if the current board fulfills the column hints.
        Only col descriptors are checked, since the rows used all
        fulfill the descriptor.
        :return: True, if the board fulfills the column hints
        """
        descriptor = self.get_col_descriptors()
        equal = descriptor == self.col_hints
        return equal

    def __find_solution(self, rowid):
        """
        Find a solution for the nonogram up to given row number.
        It tries adding all possible row configurations for the current row.
        For each iteration, if it is a valid board sofar, the next row is added (recursive)
        :param rowid: Row number of the current row
        :return: True, if a solution is found
        """
        if rowid < self.no_rows:
            for i in range(0, len(self.possible_rows[rowid])):
                self.board[rowid] = self.possible_rows[rowid][i]
                if self.match_cols_sofar(rowid + 1):
                    found = self.__find_solution(rowid + 1)
                    if found:
                        return True
        else:
            # Only check columns, only valid rows are iterated
            if self.check_board_vs_col_constraints():
                return True
        return False

    def solve(self):
        """
        Solve the given Nonogram.
        Returns a copy of the resulting board.
        The result is also stored in the internal board.
        If multiple solutions are found, only the last found solution is stored.
        :return: The solved board
        """
        self.generate_possible_rows()
        # Check to see if it faster the flip the board so we can find the solution faster
        # Flipping is determined on the number of possible iterations of the first three (and
        # last three) rows.
        if len(self.possible_rows[0])*len(self.possible_rows[1])*len(self.possible_rows[2]) > \
           len(self.possible_rows[-1])*len(self.possible_rows[-2])*len(self.possible_rows[-3]):
            self.flip()
        self.__find_solution(0)
        if self.flipped:
            self.flip()
        return self.board.copy()

    def timed_solve(self):
        """
        Solve the given Nonogram and show how long it took to solve.
        Returns a copy of the resulting board and the time taken
        The result is also stored in the internal board.
        If multiple solutions are found, only the last found solution is stored.
        :return: The solved board, time taken (ms)
        """
        start_time = time.time()
        result = self.solve()
        time_taken = time.time() - start_time
        print("--- %5.2f seconds ---" % time_taken)
        return result, time_taken

    def partial_match(self, base_list, partial):
        """
        Check if list partial is a partial match of list base_list.
        List b is a partial match if:
        - List partial is empty
        - List partial is smaller and the last element in the list  is smaller or equal to
          the same element in base_list and all other elements are equal
        :param base_list: the list to match against
        :param partial: the list to match
        :return: True, if b is a partial match of a
        """
        if len(partial) == 0:
            return True
        if len(partial) > len(base_list):
            return False
        if partial[len(partial) - 1] > base_list[len(partial) - 1]:
            return False
        for i in range(0, len(partial) - 1):
            if base_list[i] != partial[i]:
                return False
        return True

    def match_cols_sofar(self, no_rows):
        """
        Match the first rows of the board against the hints.
        It is a match if the descriptor of the first rows is a partial match
        of the hints (see function partial_match).
        :param no_rows: Number of rows to match
        :return: True, if the board so far fulfills the hints
        """
        res = []
        for colid in range(0, self.no_cols):
            res.append(self.get_col_descriptor(colid, no_rows))
        for colid in range(len(res)):
            if not self.partial_match(self.col_hints[colid], res[colid]):
                return False
        return True

    def flip(self):
        """
        Flip the board and the hints so that the last row becomes the first.
        The list of row hints is inverted
        For each column the list is inverted
        The board and solution are reverserd
        Flipping it twice results in the same set of hints
        :return:
        """
        self.flipped = not self.flipped
        self.row_hints.reverse()
        for i in self.col_hints:
            i.reverse()
        if self.possible_rows:
            self.possible_rows.reverse()
        if self.board:
            self.board.reverse()


if __name__ == '__main__':
    n = Nonogram()
    n.initialize_from_file('cactus.txt', verbose=True)
    n.timed_solve()
    n.print(fields=False, print_board=True)

    n = Nonogram()
    n.initialize_from_file('mill.txt', verbose=True)
    n.timed_solve()
    n.print(fields=False, print_board=True)

    n = Nonogram()
    n.initialize_from_file('sorbet.txt', verbose=True)
    n.timed_solve()
    n.print(fields=False, print_board=True)

    n = Nonogram()
    n.initialize_from_file('tree.txt', verbose=True)
    n.timed_solve()
    n.print(fields=False, print_board=True)

    n = Nonogram()
    n.initialize_from_file('10X6.txt', verbose=True)
    n.timed_solve()
    n.print(fields=False, print_board=True)

    n = Nonogram()
    n.initialize_from_file('15x15.txt', verbose=True)
    n.timed_solve()
    n.print(fields=False, print_board=True)
