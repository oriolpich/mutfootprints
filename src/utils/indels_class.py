from bgreference import refseq
import numpy as np


# this will classify the Indels according to the PCAWG classification
def IndelsClassifier(row, genome='hg19'):
    """
    This will classify indel and dbs events into the PCAWG categories to feed into the extraction.
    We are classifying indels based on the position in the genome we belive they have occured.

    a) In the case of deletions:
    1. if the first letter in ref is the same as in alt, we conclude that the fragment excised is REF[1:].
    2. if the first letter in ref differs to the one in alt, we conclude the entire REF has been excised.
      This means that when checking sequences, POS should not be included as in 1) since it is also deleted.

    b) In the case of insertions:
    1.if the first letter in ref is the same as in alt, we conclude that the fragment inserted is ALT[1:]
    2.if the first letter in ref differs to the one in alt, we conclude the insertion is the entire REF.

    """

    dipyr = ('C', 'T')
    comp = {'A': 'T', 'G': 'C'}
    complementary = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    dbs_list = {
        'AC_CA', 'AC_CG', 'AC_CT', 'AC_GA', 'AC_GG', 'AC_GT', 'AC_TA', 'AC_TG', 'AC_TT', 'AT_CA', 'AT_CC', 'AT_CG',
        'AT_GA', 'AT_GC', 'AT_TA', 'CC_AA', 'CC_AG', 'CC_AT', 'CC_GA', 'CC_GG', 'CC_GT', 'CC_TA', 'CC_TG', 'CC_TT',
        'CG_AT', 'CG_GC', 'CG_GT', 'CG_TA', 'CG_TC', 'CG_TT', 'CT_AA', 'CT_AC', 'CT_AG', 'CT_GA', 'CT_GC', 'CT_GG',
        'CT_TA', 'CT_TC', 'CT_TG', 'GC_AA', 'GC_AG', 'GC_AT', 'GC_CA', 'GC_CG', 'GC_TA', 'TA_AT', 'TA_CG', 'TA_CT',
        'TA_GC', 'TA_GG', 'TA_GT', 'TC_AA', 'TC_AG', 'TC_AT', 'TC_CA', 'TC_CG', 'TC_CT', 'TC_GA', 'TC_GG', 'TC_GT',
        'TG_AA', 'TG_AC', 'TG_AT', 'TG_CA', 'TG_CC', 'TG_CT', 'TG_GA', 'TG_GC', 'TG_GT', 'TT_AA', 'TT_AC', 'TT_AG',
        'TT_CA', 'TT_CC', 'TT_CG', 'TT_GA', 'TT_GC', 'TT_GG'
    }

    # ===================
    # DELETION BIG GROUP
    # ===================

    # Example --> chr1 	99072 	CT 	C
    # This means you have removed the T and you keep the C
    if row['CLASS'] == 'DEL':

        first_nuc_REF = row['REF'][0]
        first_nuc_ALT = row['ALT'][0]

        # this means it is constant
        if first_nuc_REF == first_nuc_ALT:
            flag_full = True
            endpos = len(row['REF'])
            size_del = len(row['REF'][1:])
            affected_site = row['REF'][1]
            affected_seq = row['REF'][1:]
            pos = int(row['POS'])

        else:
            size_del = len(row['REF'])
            affected_site = row['REF'][0]
            pos = int(row['POS']) - 1
            affected_seq = row['REF']
            endpos = len(row['REF']) + 1

        # FIRST CLASS : one base DELETION
        if size_del == 1:

            # if in dipyrimidines group
            left_sequence = refseq(genome, row['CHROM'], pos - 4, 5)
            right_sequence = refseq(genome, row['CHROM'], pos + 2, 5)

            # now we check for the repeats
            # we will never have count 0!
            count_eq = 1
            for i in right_sequence:
                if i == affected_site:
                    count_eq += 1
                else:
                    break
            for i in left_sequence[::-1]:
                if i == affected_site:
                    count_eq += 1
                else:
                    break
            # if we have 5 repetitions, we have this class
            if count_eq > 5:
                class_in = 'DEL_{}_1_6+'.format(comp.get(affected_site, affected_site))

            # else, we specifie the dipyr and how many repeats we have
            else:
                class_in = 'DEL_{}_1_{}'.format(comp.get(affected_site, affected_site), count_eq, )

        # SECOND CLASS: more than one base deletion
        elif size_del > 1:

            # len_ref-1 because the chunck excised also contains the nucleotide before
            seq1 = refseq(genome, row['CHROM'], pos + endpos, size_del * 5)

            # we want 3' of the reverse
            seq2 = refseq(genome, row['CHROM'], pos - size_del * 5, size_del * 5 + 1)[::-1]

            count_eq = 1

            # split the sequence into bins of the same size of the deleted region
            # Right location
            splitted = [seq1[i:i + size_del] for i in range(0, len(seq1), size_del)]
            for i in splitted:
                if i == affected_seq:
                    count_eq += 1
                else:
                    break

            # Left location
            splitted = [seq2[i:i + size_del] for i in range(0, len(seq2), size_del)]
            for i in splitted:
                if i == affected_seq[::-1]:
                    count_eq += 1
                else:
                    break

            # if the count is equal or greater than 5, we have this class
            if count_eq > 5:

                # this yields the DEL_repeats_5+_5+
                if size_del >= 5:
                    class_in = 'DEL_repeats_5+_6+'
                else:
                    class_in = 'DEL_repeats_{}_6+'.format(size_del)
            else:

                # if we have some repeats found and they are less than 5,
                # then they belong to the next class
                if count_eq > 1:
                    if size_del >= 5:
                        class_in = 'DEL_repeats_5+_{}'.format(count_eq)
                    else:
                        class_in = 'DEL_repeats_{}_{}'.format(size_del, count_eq)

                # if no full repeat is found, then give opportunity to microhomology
                else:

                    # get sequence on the right
                    right_seq = refseq(genome, row['CHROM'], pos + endpos, size_del)
                    # get sequence on the left
                    left_seq = refseq(genome, row['CHROM'], pos - size_del, size_del + 1)

                    good = 0

                    # we go down the size of the indel
                    for i in np.arange(size_del - 1, 0, -1):

                        # check right side
                        tocheck = affected_seq[:i]
                        tocheck_right = right_seq[:i]
                        if tocheck == tocheck_right:
                            good = i
                            break

                        # check the left side
                        tocheck = affected_seq[::-1][:i]
                        tocheck_left = left_seq[::-1][:i]

                        if tocheck == tocheck_left:
                            good = i
                            break

                    # if microhomology has been detected
                    if good > 0:
                        # this yields the DEL_repeats_5+_5+
                        if good >= 5:
                            good = '5+'

                        if size_del >= 5:
                            class_in = 'DEL_MH_5+_{}'.format(good)
                        else:
                            class_in = 'DEL_MH_{}_{}'.format(size_del, good)

                    # else this means a deletion with 0 repetitions
                    else:
                        if size_del >= 5:
                            class_in = 'DEL_repeats_5+_1'  # we put one according to their info
                        else:
                            class_in = 'DEL_repeats_{}_1'.format(size_del)

    # ===================
    # INSERTIONS BIG GROUP
    # ===================
    elif row['CLASS'] == 'INS':

        first_nuc_REF = row['REF'][0]
        first_nuc_ALT = row['ALT'][0]

        if first_nuc_REF == first_nuc_ALT:
            flag_full = True
            endpos = len(row['REF'])
            size_del = len(row['ALT'][1:])
            affected_site = row['ALT'][1]
            affected_seq = row['ALT'][1:]
            pos = int(row['POS'])

        else:
            flag_full = False
            size_del = len(row['ALT'])
            affected_site = row['ALT'][0]
            pos = int(row['POS']) - 1
            affected_seq = row['ALT']
            endpos = len(row['REF']) + 1

        # FIRST CLASS : one base INSERTION
        if size_del == 1:

            # if in dipyrimidines group
            # pos+2 because the deletion is mapped just at the beginning of POS!
            right_sequence = refseq(genome, row['CHROM'], pos + endpos, 5)

            # else we want whathever is on the 3' (because we would reverse it)
            # we want 3' of the reversed. This should include the first nucleotide!
            left_sequence = refseq(genome, row['CHROM'], pos - 4, 5)

            # we will never have count 0!
            count_eq = 0

            for i in right_sequence:
                if i == affected_site:
                    count_eq += 1
                else:
                    break
            for i in left_sequence[::-1]:
                if i == affected_site:
                    count_eq += 1
                else:
                    break

            if count_eq >= 5:
                class_in = 'INS_{}_1_5+'.format(comp.get(affected_site, affected_site))

            # else, we specifie the dipyr and how many repeats we have
            else:
                class_in = 'INS_{}_1_{}'.format(comp.get(affected_site, affected_site), count_eq, )

        elif size_del > 1:

            # len_ref-1 because the chunck excised also contains the nucleotide before
            seq1 = refseq(genome, row['CHROM'], pos + endpos, size_del * 5)
            seq2 = refseq(genome, row['CHROM'], pos - size_del * 5, size_del * 5 + 1)[::-1]

            count_eq = 0

            # split the sequence into bins of the same size of the deleted region
            splitted = [seq1[i:i + size_del] for i in range(0, len(seq1), size_del)]

            for i in splitted:
                if i == affected_seq:
                    count_eq += 1
                else:
                    break

            splitted = [seq2[i:i + size_del] for i in range(0, len(seq2), size_del)]
            for i in splitted:
                if i == affected_seq[::-1]:
                    count_eq += 1
                else:
                    break
            # if the count is equal or greater than 5, we have this class
            if count_eq >= 5:

                # this yields the INS_repeats_5+_5+
                if size_del >= 5:
                    class_in = 'INS_repeats_5+_5+'
                else:
                    class_in = 'INS_repeats_{}_5+'.format(size_del)

            else:
                if size_del >= 5:
                    class_in = 'INS_repeats_5+_{}'.format(count_eq)
                else:
                    class_in = 'INS_repeats_{}_{}'.format(size_del, count_eq)

    # ===================
    # DBS BIG GROUP
    # ===================
    elif row['CLASS'] == 'DBS':

        # When merging the reverse complementary doublet base substitutions classes into one class, 12 of mutation
        # classes have no strandness (e.g. CG > AT), resulting in 78 classes of doublet base substitutions

        class_in = '{}_{}'.format(row['REF'], row['ALT'])
        if class_in not in dbs_list:

            # AC>CA	GT>TG
            class_in = '{}{}_{}{}'.format(
                complementary[row['REF'][1]], complementary[row['REF'][0]],
                complementary[row['ALT'][1]], complementary[row['ALT'][0]]
            )

    else:
        class_in = 'NOTGOOD'

    return class_in
