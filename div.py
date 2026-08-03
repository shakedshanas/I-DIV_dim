from numpy import array
import math

def DIV(arr, spreading_percentage, tpm_or_count, ceil=False, min_eve=0):
    """
    Adaptation of Shannon's eveness index.
    Computes array's isoform diversity.
    
    Args:
        arr: An array of isoform expression.
        spreading_percentage: An int which decides the spreading percentage we want to spread.
        tpm_or_count: A string which decides whether to read TPM files or Transcript Read Counts file.
        ceil: A boolean which decides whether to use ceiling.
        min_eve: A boolean which decides whether to calculate the minimum diversity.
    
    Returns:
        An number which represent the diversity of isoforms in a selected gene.
    """
    
    # Calculating how many actual splice sites in total, meaning sites with readings > 0. We want that len will be > 1
    # We don't want to use spreading on genes with 1 isoforms, this will return a value with maximum diversity - we return 0 diversity instead.
    if len(array(arr)[array(arr)>0]) <= 1:
        return 0   
    
    # Calculating the percentage of read counts we want to spread 
    p = ((sum(arr)*(spreading_percentage/100)) / len(arr))
    
    # Decides whether to use ceiling or not.
    # If we:
    #   use count data, 
    #   define ceil==true 
    #   and 0<p<1 (greater than 0 in case we want no spreading), 
    # then:
    #   ceiling is applied
    if tpm_or_count == 'count' and ceil and (0<p<1):
        p = 1
    

    if min_eve: # If true, then we calcualte the minimum diversity
        ''' minimum diversity: one isoform that contains all the reads in it. 
        for example: [1000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        we take the minimum diversity of an arary. 
        for example: [250, 250, 250, 250, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        and we measure the minimum diversity as it it has one isoform that contains all the reads in the array
        like in the first example. The first eleemnt will have the sum of reads.
        we add the p (spreading percentage) let's say that p=1%
        for example: [1000+10, 0+10, 0+10, 0+10, 0+10, 0+10, 0+10, 0+10, 0+10, 0+10, 0+10, 0+10, 0+10]
        '''
        arr_new = array([sum(arr)+p] + [p for _ in range(1,len(arr))] )
    else:
        arr_new = array(arr) + p # spreading
    
    
    possible_ss = len(arr_new) #computing how many POSSIBLE splice sites (including sites with 0 readings)
    total_reads = sum(arr_new) #computing how many readings in TOTAL.

    
    #computing shannon's eveness index
    index_value = 0
    for splice_site in arr_new:
        x = splice_site/(total_reads)
        if x!=0:
            index_value += (x)*(math.log(x)) # default base is e
   
    return (round(-1*index_value, 3) / math.log(possible_ss)) # default base is e

def DIV_min(arr, spreading_percentage, tpm_or_count, ceil=False):
    """
    Computes array's (isoforms splicing choices in a gene) MINIMUM diversity.
    
    It is important to note that we DON'T calculate the real minimum (obviously it would be 0).
    All the elements in the array are summed, and a new array is defined in which the first element contains
    the sum of all elements, and all other elements contains ONLY their portion of the spreading percentage
    """
    
    return DIV(arr, spreading_percentage, tpm_or_count, ceil, min_eve=1)


def filterrows(rows_pd, column_filter_by,gid):
    ''' GTEx matrices are too big to load all at once, so we filter the rows by gid.
    This function returns the rows that should be skipped when loading the data.
    
    Args:
        rows_pd: row numbers of dataframe.
        column_filter_by: The column name to filter by (column name of gene IDs).
        gid: The gene ID to filter for.
    
    Returns:
        A list of row indices to skip when loading the data.
    '''
    # rows_pd_stay = rows_pd[rows_pd[column_filter_by].str.contains(gid)].index
    rows_pd_stay = rows_pd.loc[rows_pd[column_filter_by] == gid].index
    # print(f"Rows to stay: {rows_pd_stay}")
    rows_pd_skip = [x + 1  for x in rows_pd.index if x not in rows_pd_stay] 
    
    return rows_pd_skip