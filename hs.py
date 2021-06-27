from hsCommand import *

def main():
    chroms_name = []
    indexes_chrom_offset = []
    pixels_bin1_id = []
    pixels_bin2_id = []
    pixels_count = []

    actual_contigs = []

    while True:
        command = input().split()

        if command[0] == "read":
            chroms_name, indexes_chrom_offset, pixels_bin1_id, pixels_bin2_id, pixels_count = read_hic(command[1])
            
            for i in range(len(indexes_chrom_offset) - 1):
                actual_contigs.append({"contig": i, "ori_start": indexes_chrom_offset[i], 
                                       "ori_end": indexes_chrom_offset[i + 1], "curr_start": indexes_chrom_offset[i], 
                                       "curr_end": indexes_chrom_offset[i + 1], 
                                       "length": indexes_chrom_offset[i + 1] - indexes_chrom_offset[i],
                                       "orientation": "forward"})

        elif command[0] == "move":
            move_contig(actual_contigs, int(command[1]), int(command[2]))

        elif command[0] == "reverse":
            reverse_contig(actual_contigs, int(command[1]))

        elif command[0] == "savepng":
            save_png(command[1], int(command[2]), int(command[3]), int(command[4]), int(command[5]), actual_contigs, pixels_bin1_id, pixels_bin2_id, pixels_count)

        elif command[0] == "savefasta":
            save_fasta(command[1], command[2], actual_contigs, chroms_name, )

        elif command[0] == "exit":
            exit()

        elif command[0] == "help":
            print("Usage:",
                  "read <path-to-cool-file>",
                  "move <num-of-contig> <where>",
                  "reverse <num-of-contig>",
                  "savepng <name-of-the-file> <i1> <i2> <j1> <j2>",
                  "savefasta <path-to-original-fasta> <edited-fasta>", sep = "\n")

        else:
            print("No such command")
            


if __name__ == "__main__":
    main()