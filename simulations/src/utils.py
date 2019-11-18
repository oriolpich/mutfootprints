""" tri-nucleotide encoding conversion utils """


def key_func1(r):

    subs = r['Type']
    triplet = r['SubType']

    return subs[0] + subs[-1], triplet[0] + triplet[-1]


def key_func2(r):

    subs, triplet = tuple(r['feature'].split('_at_'))

    return subs[0] + subs[-1], triplet[0] + triplet[-1]
