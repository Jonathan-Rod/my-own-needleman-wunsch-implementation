import sys
import csv

# constants given in instructions
MATCH, MISMATCH, GAP = 1, -1, -2


def get_backtracking(seq_1: str, seq_2: str, btrack_grid: list) -> str:
    """Generate the global alignment sequences

    Arguments:
        seq_1 (str): first sequence
        seq_2 (str): second sequence
        btrack_grid(list): character matrix

    Returns:
        str: alignment sequences
    """
    alig_1, alig_2 = '', ''
    i, j = len(btrack_grid) - 1, len(btrack_grid[0]) - 1
    char = btrack_grid[i][j]
    while char != None:
        # diagonal case:
        #   add the character in two alignment sequences
        #   moves diagonally backward
        if char == 'D':
            alig_1 = seq_1[i - 1] + alig_1
            alig_2 = seq_2[j - 1] + alig_2
            i -= 1
            j -= 1
        # left case:
        #   add the gap in first alignment sequence
        #   add the character in second alignment sequence
        #   moves left
        elif char == 'L':
            alig_1 = '-' + alig_1
            alig_2 = seq_2[j - 1] + alig_2
            j -= 1
        # up case:
        #   add the character in first alignment sequence
        #   add the gap in second alignment sequence
        #   moves up
        elif char == 'U':
            alig_1 = seq_1[i - 1] + alig_1
            alig_2 = '-' + alig_2
            i -= 1
        char = btrack_grid[i][j]
    return alig_1 + ' ' + alig_2


def compute(seq_1: str, seq_2: str, grid: list, btrack_grid: list):
    """Fill's the scoring matrix using the Needleman Wunsh algorithm.

    Arguments:
        seq_1 (str): first sequence
        seq_2 (str): second sequence
        grid(list): scoring matrix
        btrack_grid(list): character matrix
    """
    for i in range(1, len(grid)):
        char_1 = seq_1[i-1:i]
        for j in range(1, len(grid[i])):
            char_2 = seq_2[j-1:j]
            # Based on characters, gap, and values on scoring matrix
            choice_1 = grid[i-1][j-1] + \
                (MATCH if char_1 == char_2 else MISMATCH)
            choice_2 = grid[i-1][j] + GAP
            choice_3 = grid[i][j-1] + GAP
            # Set the maximum value obtained
            maximum = grid[i][j] = max(choice_1, choice_2, choice_3)
            char = None
            if maximum == choice_3:
                char = 'L'
            elif maximum == choice_2:
                char = 'U'
            elif maximum == choice_1:
                char = 'D'
            # Set the character based on maximum value (left, diagonal, or up)
            btrack_grid[i][j] = char


def set_initialization_step(grid: list, btrack_grid: list):
    """Initialize scoring matrix and direction matrix

    Arguments:
        grid(list): scoring matrix
        btrack_grid(list): character matrix
    """
    for j in range(len(grid[0])):
        grid[0][j] = GAP * j
        btrack_grid[0][j] = 'L'
    for i in range(len(grid)):
        grid[i][0] = GAP * i
        btrack_grid[i][0] = 'U'
    btrack_grid[0][0] = None


def get_grid(rows: int, cols: int) -> list:
    """Creates an empty array with a specified number of rows and columns

    Arguments:
        rows (int): number of rows
        cols (int): number of columns

    Returns:
        list: empty array
    """

    new_grid = []
    for _ in range(rows):
        temp_row = []
        for _ in range(cols):
            temp_row.append(None)
        new_grid.append(temp_row)
    return new_grid


def get_outputs(inputs: list) -> list:
    """Returns a list of alignment sequences

    arguments:
        inputs(list): sequences to compute

    Returns:
        list: alignment sequences
    """
    outputs = []
    for row in inputs:
        seq_1, seq_2 = row[0], row[1]
        row_size, col_size = len(seq_1) + 1, len(seq_2) + 1
        grid = get_grid(row_size, col_size)
        btrack_grid = get_grid(row_size, col_size)
        set_initialization_step(grid, btrack_grid)
        compute(seq_1, seq_2, grid, btrack_grid)
        seq_3 = get_backtracking(seq_1, seq_2, btrack_grid)
        score = str(grid[row_size - 1][col_size - 1])
        outputs.append(seq_3 + ' ' + score)
    return outputs


def get_csv_inputs(csv_path: str) -> list:
    """Returns a list with the sequences in the csv file

    Arguments:
        csv_path(str): csv file

    Returns:
        list: sequences
    """
    rows = []
    with open(csv_path, 'r') as file:
        csvreader = list(csv.reader(file))
        for i in range(1, len(csvreader)):
            rows.append(csvreader[i])
    return rows


def main():
    """Needleman-Wunsch implementation
    """
    # extracts path from sys argv
    path = sys.argv[len(sys.argv) - 1]
    # extracts adn sequences from input csv file
    inputs = get_csv_inputs(path)
    # standard output line
    outputs = get_outputs(inputs)
    for row in outputs:
        print(row)


if len(sys.argv) > 1:
    main()
