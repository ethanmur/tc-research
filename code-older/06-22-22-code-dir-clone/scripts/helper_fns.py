import os

def display_data_files( path, data_type, print_status):
    """
    This function simply prints out and returns a list of data files found at a
    specified path.
    :param path: the user specified path to data files. This argument must always
                 be present.
    :param data_type: must include a data identifier, like 'crl', 'tdr', 'in-situ',
                  or 'dropsonde'. Can also include 'hide-list' to surpress printing
                  a list of the data
    :param print_status: Implicitly print the file names if given 'show-list',
                  unless "hide-list" is passed instead.
    :return: A list of file names found within the specified folder.
    """

    file_names = []
    for (dirpath, dirnames, file) in os.walk( path):
        file_names.extend(file)
        break

    if data_type == 'crl' and print_status == 'show-list':
        print( "crl data files:")
        for number in range( len( file_names)):
            print( str( number) + ") " + file_names[ number])
    elif data_type == 'tdr' and print_status == 'show-list':
        print( "tdr data files:")
        for number in range( len( file_names)):
            print( str( number) + ") " + file_names[ number])
    elif data_type == 'in-situ' and print_status == 'show-list':
        print( "in situ data files:")
        for number in range( len( file_names)):
            print( str( number) + ") " + file_names[ number])
    elif data_type == 'dropsonde' and print_status == 'show-list':
        print( "dropsonde data files:")
        for number in range( len( file_names)):
            print( str( number) + ") " + file_names[ number])
    elif data_type == 'goes' and print_status == 'show-list':
        print( "GOES satellite data files:")
        for number in range( len( file_names)):
            print( str( number) + ") " + file_names[ number])

    return file_names
