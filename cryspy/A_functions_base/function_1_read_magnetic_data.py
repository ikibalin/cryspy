# ORIGINAL FORTRAN PROGRAM
#       program read_magnetic_data
# * read data about magnetic space groups
# * input data from magnetic_table.dat
#       implicit none
#       integer i,j,k,m,n

# * for the ith nonhexagonal point operator:
# * point_op_label(i): point operator symbol (from Litvin)
#       character point_op_label(48)*8
# * point_op_xyz(i): point operator in x,y,z notation
#       character point_op_xyz(48)*10
# * point_op_matrix(i): point operator matrix
#       integer point_op_matrix(3,3,48)

# * for the ith hexagonal point operator:
# * point_op_hex_label(i): point operator symbol (from Litvin)
#       character point_op_hex_label(24)*8
# * point_op_hex_xyz(i): point operator in x,y,z notation
#       character point_op_hex_xyz(24)*10
# * point_op_hex_matrix(i): point operator matrix
#       integer point_op_hex_matrix(3,3,24)

# * number of magnetic space groups
#       integer magcount
#       parameter(magcount=1651)

# * for the ith magnetic space group
# * nlabel_bns(i): numerical label in BNS setting
#       character nlabel_bns(magcount)*12
# * nlabel_parts_bns(j,i): jth part of nlabel_bns
#       integer nlabelparts_bns(2,magcount)
# * label_bns(i): group symbol
#       character spacegroup_label_bns(magcount)*14
# * nlabel_og(i): numerical label in OG setting
#       character nlabel_og(magcount)*12
# * nlabel_parts_og(j,i): jth part of nlabel_og
#       integer nlabelparts_og(3,magcount)
# * label_og(i): group symbol
#       character spacegroup_label_og(magcount)*14
# * magtype(i): type of magnetic space group (1-4)
#       integer magtype(magcount)

# * BNS-OG transformation (if type-4)
# * bnsog_point_op(j,k,i): 3x3 point operator part of transformation
#       integer bnsog_point_op(3,3,magcount)
# * bnsog_origin(j,i): translation part of transformation
# * bnsog_point_origin(i): common denominator
#       integer bnsog_origin(3,magcount)
#       integer bnsog_origin_denom(magcount)

# * iops_count(i): number of point operators
#       integer ops_count(magcount)
# * wyckoff_count(i): number of wyckoff sites
#       integer wyckoff_site_count(magcount)
# * wyckoff_pos_count(j,i): number of positions in jth wyckoff site
#       integer wyckoff_pos_count(27,magcount)
# * wyckoff_mult(j,i): multiplicity for jth wyckoff site
#       integer wyckoff_mult(27,magcount)
# * wyckoff_label(j,i): symbol (a,b,c,...,z,alpha) for jth wyckoff site
#       character wyckoff_label(27,magcount)

# * for BNS setting
# * lattice_bns_vectors_count(i): number of lattice vectors defining the lattice
#       integer lattice_bns_vectors_count(magcount)
# * lattice_bns_vectors(k,j,i): kth component of the jth lattice vector
# * lattice_bns_vectors_denom(j,i): common denominator
#       integer lattice_bns_vectors(3,6,magcount)
#       integer lattice_bns_vectors_denom(6,magcount)
# * for jth operator
# * ops_bns_point_op(j,i): point operator part
#       integer ops_bns_point_op(96,magcount)
# * ops_bns_trans(k,j,i): kth component of translation part
# * ops_bns_trans_denom(j,i): common denominator
#       integer ops_bns_trans(3,96,magcount)
#       integer ops_bns_trans_denom(96,magcount)
# * ops_bns_timeinv(j,i): 1=no time inversion, -1=time inversion
#       integer ops_bns_timeinv(96,magcount)
# * for jth wyckoff site
# * wyckoff_bns_fract(k,j,i): kth component of fractional part of wyckoff position
# * wyckoff_bns_fract_denom(j,i): common denominator
#       integer wyckoff_bns_fract(3,96,27,magcount)
#       integer wyckoff_bns_fract_denom(96,27,magcount)
# * wyckoff_bns_xyz(m,k,j,i): mth component to coeffcient of kth paramater (x,y,z)
#       integer wyckoff_bns_xyz(3,3,96,27,magcount)
# * wyckoff_bns_mag(m,k,j,i): mth component to coeffcient of kth magnetic
# * paramater (mx,my,mz)
#       integer wyckoff_bns_mag(3,3,96,27,magcount)

# * for OG setting (for type-4 groups)
# * lattice_og_vectors_count(i): number of lattice vectors defining the lattice
#       integer lattice_og_vectors_count(magcount)
# * lattice_og_vectors(k,j,i): kth component of the jth lattice vector
# * lattice_og_vectors_denom(j,i): common denominator
#       integer lattice_og_vectors(3,6,magcount)
#       integer lattice_og_vectors_denom(6,magcount)
# * for jth operator
# * ops_og_point_op(j,i): point operator part
#       integer ops_og_point_op(96,magcount)
# * ops_og_trans(k,j,i): kth component of translation part
# * ops_og_trans_denom(j,i): common denominator
#       integer ops_og_trans(3,96,magcount)
#       integer ops_og_trans_denom(96,magcount)
# * ops_og_timeinv(j,i): 1=no time inversion, -1=time inversion
#       integer ops_og_timeinv(96,magcount)
# * for jth wyckoff site
# * wyckoff_og_fract(k,j,i): kth component of fractional part of wyckoff position
# * wyckoff_og_fract_denom(j,i): common denominator
#       integer wyckoff_og_fract(3,96,27,magcount)
#       integer wyckoff_og_fract_denom(96,27,magcount)
# * wyckoff_og_xyz(m,k,j,i): mth component to coeffcient of kth paramater (x,y,z)
#       integer wyckoff_og_xyz(3,3,96,27,magcount)
# * wyckoff_og_mag(m,k,j,i): mth component to coeffcient of kth magnetic
# * paramater (mx,my,mz)
#       integer wyckoff_og_mag(3,3,96,27,magcount)

# ******************************************************************************
# * open data file
#       open(30,file='magnetic_data.txt')
# * read nonhexangonal point operators
#       do i=1,48
#         read(30,*)n,point_op_label(i),point_op_xyz(i),
#      $       ((point_op_matrix(k,j,i),j=1,3),k=1,3)
#         if(n.ne.i)stop
#      $       'error in numbering of nonhexagonal point operators'
#       enddo
# * read hexangonal point operators
#       do i=1,24
#         read(30,*)n,point_op_hex_label(i),
#      $       point_op_hex_xyz(i),
#      $       ((point_op_hex_matrix(k,j,i),j=1,3),k=1,3)
#         if(n.ne.i)stop
#      $       'error in numbering of hexagonal point operators'
#       enddo
# * read data for each magnetic space group
#       do i=1,1651
#         read(30,*)(nlabelparts_bns(j,i),j=1,2),nlabel_bns(i),
#      $       spacegroup_label_bns(i),(nlabelparts_og(j,i),j=1,3),
#      $       nlabel_og(i),spacegroup_label_og(i)
#         read(30,*)magtype(i)
#         if(magtype(i).eq.4)then
#           read(30,*)((bnsog_point_op(j,k,i),j=1,3),k=1,3),
#      $         (bnsog_origin(j,i),j=1,3),bnsog_origin_denom(i)
#         endif
#         read(30,*)ops_count(i)
#         read(30,*)(ops_bns_point_op(j,i),(ops_bns_trans(k,j,i),k=1,3),
#      $       ops_bns_trans_denom(j,i),ops_bns_timeinv(j,i),
#      $       j=1,ops_count(i))
#         read(30,*)lattice_bns_vectors_count(i)
#         read(30,*)((lattice_bns_vectors(k,j,i),k=1,3),
#      $       lattice_bns_vectors_denom(j,i),
#      $       j=1,lattice_bns_vectors_count(i))
#         read(30,*)wyckoff_site_count(i)
#         do j=1,wyckoff_site_count(i)
#           read(30,*)wyckoff_pos_count(j,i),wyckoff_mult(j,i),
#      $         wyckoff_label(j,i)
#           do k=1,wyckoff_pos_count(j,i)
#             read(30,*)(wyckoff_bns_fract(m,k,j,i),m=1,3),
#      $           wyckoff_bns_fract_denom(k,j,i),
#      $           ((wyckoff_bns_xyz(m,n,k,j,i),m=1,3),n=1,3),
#      $           ((wyckoff_bns_mag(m,n,k,j,i),m=1,3),n=1,3)
#           enddo
#         enddo
#         if(magtype(i).eq.4)then
#         read(30,*)ops_count(i)
#         read(30,*)(ops_og_point_op(j,i),(ops_og_trans(k,j,i),k=1,3),
#      $       ops_og_trans_denom(j,i),ops_og_timeinv(j,i),
#      $       j=1,ops_count(i))
#         read(30,*)lattice_og_vectors_count(i)
#         read(30,*)((lattice_og_vectors(k,j,i),k=1,3),
#      $       lattice_og_vectors_denom(j,i),
#      $       j=1,lattice_og_vectors_count(i))
#         read(30,*)wyckoff_site_count(i)
#         do j=1,wyckoff_site_count(i)
#           read(30,*)wyckoff_pos_count(j,i),wyckoff_mult(j,i),
#      $         wyckoff_label(j,i)
#           do k=1,wyckoff_pos_count(j,i)
#             read(30,*)(wyckoff_og_fract(m,k,j,i),m=1,3),
#      $           wyckoff_og_fract_denom(k,j,i),
#      $           ((wyckoff_og_xyz(m,n,k,j,i),m=1,3),n=1,3),
#      $           ((wyckoff_og_mag(m,n,k,j,i),m=1,3),n=1,3)
#           enddo
#         enddo
#         endif
#       enddo
# * close data file
#       close(30)
#       end

import os
import numpy

# for the ith nonhexagonal point operator:
# point_op_label(i): point operator symbol (from Litvin)
point_op_label = numpy.zeros(shape = (48, ), dtype='<8U')

# point_op_xyz(i): point operator in x,y,z notation
point_op_xyz = numpy.zeros(shape = (48, ), dtype='<10U')

# point_op_matrix(i): point operator matrix
point_op_matrix = numpy.zeros(shape = (3, 3, 48, ), dtype=int)

# for the ith hexagonal point operator:
# point_op_hex_label(i): point operator symbol (from Litvin)
point_op_hex_label = numpy.zeros(shape = (24, ), dtype='<8U')

# point_op_hex_xyz(i): point operator in x,y,z notation
point_op_hex_xyz = numpy.zeros(shape = (24, ), dtype='<10U')

# point_op_hex_matrix(i): point operator matrix
point_op_hex_matrix = numpy.zeros(shape = (3, 3, 24, ), dtype=int)

      
# number of magnetic space groups
magcount = 1651

# for the ith magnetic space group
# nlabel_bns(i): numerical label in BNS setting
nlabel_bns = numpy.zeros(shape = (magcount, ), dtype='<12U')

# nlabel_parts_bns(j,i): jth part of nlabel_bns
nlabelparts_bns = numpy.zeros(shape = (2, magcount, ), dtype=int)

# label_bns(i): group symbol
spacegroup_label_bns =  numpy.zeros(shape = (magcount, ), dtype='<14U')

# nlabel_og(i): numerical label in OG setting
nlabel_og =  numpy.zeros(shape = (magcount, ), dtype='<12U')

# nlabel_parts_og(j,i): jth part of nlabel_og
nlabelparts_og = numpy.zeros(shape = (3, magcount, ), dtype=int)

# label_og(i): group symbol
spacegroup_label_og =  numpy.zeros(shape = (magcount, ), dtype='<14U')

# magtype(i): type of magnetic space group (1-4)
magtype = numpy.zeros(shape = (magcount, ), dtype=int)

# BNS-OG transformation (if type-4)
# bnsog_point_op(j,k,i): 3x3 point operator part of transformation
bnsog_point_op = numpy.zeros(shape = (3, 3, magcount, ), dtype=int)

# bnsog_origin(j,i): translation part of transformation
# bnsog_point_origin(i): common denominator
bnsog_origin = numpy.zeros(shape = (3, magcount, ), dtype=int)
bnsog_origin_denom = numpy.zeros(shape = (magcount, ), dtype=int)

# iops_count(i): number of point operators
ops_count = numpy.zeros(shape = (magcount, ), dtype=int)

# wyckoff_count(i): number of wyckoff sites
wyckoff_site_count = numpy.zeros(shape = (magcount, ), dtype=int)

# wyckoff_pos_count(j,i): number of positions in jth wyckoff site
wyckoff_pos_count = numpy.zeros(shape = (27, magcount, ), dtype=int)

# wyckoff_mult(j,i): multiplicity for jth wyckoff site
wyckoff_mult = numpy.zeros(shape = (27, magcount, ), dtype=int)
      
# wyckoff_label(j,i): symbol (a,b,c,...,z,alpha) for jth wyckoff site
wyckoff_label = numpy.zeros(shape = (27, magcount, ), dtype="<1U")

# for BNS setting
# lattice_bns_vectors_count(i): number of lattice vectors defining the lattice
lattice_bns_vectors_count = numpy.zeros(shape = (magcount, ), dtype=int)

# lattice_bns_vectors(k,j,i): kth component of the jth lattice vector
# lattice_bns_vectors_denom(j,i): common denominator
lattice_bns_vectors = numpy.zeros(shape = (3, 6, magcount, ), dtype=int)
lattice_bns_vectors_denom = numpy.zeros(shape = (6, magcount, ), dtype=int)

# for jth operator
# ops_bns_point_op(j,i): point operator part
ops_bns_point_op = numpy.zeros(shape = (96, magcount, ), dtype=int)

# ops_bns_trans(k,j,i): kth component of translation part
# ops_bns_trans_denom(j,i): common denominator
ops_bns_trans = numpy.zeros(shape = (3, 96, magcount, ), dtype=int)
ops_bns_trans_denom = numpy.zeros(shape = (96, magcount, ), dtype=int)

# ops_bns_timeinv(j,i): 1=no time inversion, -1=time inversion
ops_bns_timeinv = numpy.zeros(shape = (96, magcount, ), dtype=int)

# for jth wyckoff site
# wyckoff_bns_fract(k,j,i): kth component of fractional part of wyckoff position
# wyckoff_bns_fract_denom(j,i): common denominator
wyckoff_bns_fract = numpy.zeros(shape = (3, 96, 27, magcount, ), dtype=int)
wyckoff_bns_fract_denom = numpy.zeros(shape = (96, 27, magcount, ), dtype=int)

# wyckoff_bns_xyz(m,k,j,i): mth component to coeffcient of kth paramater (x,y,z)
wyckoff_bns_xyz = numpy.zeros(shape = (3, 3, 96, 27, magcount, ), dtype=int)

# wyckoff_bns_mag(m,k,j,i): mth component to coeffcient of kth magnetic
# paramater (mx,my,mz)
wyckoff_bns_mag = numpy.zeros(shape = (3, 3, 96, 27, magcount, ), dtype=int)


# for OG setting (for type-4 groups)
# lattice_og_vectors_count(i): number of lattice vectors defining the lattice
lattice_og_vectors_count = numpy.zeros(shape = (magcount, ), dtype=int)

# lattice_og_vectors(k,j,i): kth component of the jth lattice vector
# lattice_og_vectors_denom(j,i): common denominator
lattice_og_vectors = numpy.zeros(shape = (3, 6, magcount, ), dtype=int)
lattice_og_vectors_denom = numpy.zeros(shape = (6, magcount, ), dtype=int)

# for jth operator
# ops_og_point_op(j,i): point operator part
ops_og_point_op = numpy.zeros(shape = (96, magcount, ), dtype=int)

# ops_og_trans(k,j,i): kth component of translation part
# ops_og_trans_denom(j,i): common denominator
ops_og_trans = numpy.zeros(shape = (3, 96, magcount, ), dtype=int)
ops_og_trans_denom = numpy.zeros(shape = (96, magcount, ), dtype=int)

# ops_og_timeinv(j,i): 1=no time inversion, -1=time inversion
ops_og_timeinv = numpy.zeros(shape = (96, magcount, ), dtype=int)

# for jth wyckoff site
# wyckoff_og_fract(k,j,i): kth component of fractional part of wyckoff position
# wyckoff_og_fract_denom(j,i): common denominator
wyckoff_og_fract = numpy.zeros(shape = (3, 96, 27, magcount, ), dtype=int)
wyckoff_og_fract_denom = numpy.zeros(shape = (96, 27, magcount, ), dtype=int)

# wyckoff_og_xyz(m,k,j,i): mth component to coeffcient of kth paramater (x,y,z)
wyckoff_og_xyz = numpy.zeros(shape = (3, 3, 96, 27, magcount, ), dtype=int)

# wyckoff_og_mag(m,k,j,i): mth component to coeffcient of kth magnetic
# paramater (mx,my,mz)
wyckoff_og_mag = numpy.zeros(shape = (3, 3, 96, 27, magcount, ), dtype=int)


F_MAG_DATA = os.path.join(os.path.dirname(__file__), "magnetic_data.txt")
      
def read_magnetic_data():
    with open(F_MAG_DATA, "r") as fid:
        l_cont = fid.readlines()
    # read nonhexangonal point operators
    i_line = 0
    for i in range(48):
        l_h  = l_cont[i_line].split(); i_line += 1
        n = int(l_h[0])
        point_op_label[i] = l_h[1]
        point_op_xyz[i] = l_h[2]
        point_op_matrix[:,:,i] = numpy.array(l_h[3:3+9], dtype=int).reshape(
            3, 3)
        if n != (i+1):
            break
    # read hexangonal point operators
    for i in range(24):
        l_h  = l_cont[i_line].split(); i_line += 1
        n = int(l_h[0])
        point_op_hex_label[i] = l_h[1]
        point_op_hex_xyz[i] = l_h[2]
        point_op_hex_matrix[:,:,i] = numpy.array(l_h[3:3+9], dtype=int).reshape(
            3, 3)
        if n != (i+1):
            break

    for i in range(magcount):
        l_h  = l_cont[i_line].split(); i_line += 1
        nlabelparts_bns[:,i] = numpy.array(l_h[0:2], dtype=int)
        nlabel_bns[i] = l_h[2]
        spacegroup_label_bns[i] = l_h[3]
        nlabelparts_og[:, i] = numpy.array(l_h[4:7], dtype=int)
        nlabel_og[i] = l_h[7]
        spacegroup_label_og[i] = l_h[8]

        l_h  = l_cont[i_line].split(); i_line += 1
        magtype[i] = l_h[0]
        
        if magtype[i] == 4:
            l_h  = l_cont[i_line].split(); i_line += 1
            bnsog_point_op[:,:,i] = numpy.transpose(
                numpy.array(l_h[0:9], dtype=int).reshape(3, 3))
            bnsog_origin[:, i] = numpy.array(l_h[9:9+3], dtype=int)
            bnsog_origin_denom[i] = l_h[12]
            
        l_h  = l_cont[i_line].split(); i_line += 1
        ops_count[i] = l_h[0]

        for j in range(ops_count[i]):
            if j%4 == 0:
                l_h  = l_cont[i_line].split(); i_line += 1
            ops_bns_point_op[j, i], ops_bns_trans[0,j,i], \
                ops_bns_trans[1,j,i], ops_bns_trans[2,j,i], \
                ops_bns_trans_denom[j,i], ops_bns_timeinv[j,i] = \
                l_h[0+6*j%4], l_h[1+6*j%4], l_h[2+6*j%4], l_h[3+6*j%4], \
                l_h[4+6*j%4], l_h[5+6*j%4]

        l_h  = l_cont[i_line].split(); i_line += 1
        lattice_bns_vectors_count[i] = l_h[0]

        l_h  = l_cont[i_line].split(); i_line += 1
        for j in range(lattice_bns_vectors_count[i]):
            lattice_bns_vectors[0, j, i], lattice_bns_vectors[1, j, i], \
                lattice_bns_vectors[2, j, i], \
                lattice_bns_vectors_denom[j, i] = \
                l_h[0+4*j], l_h[1+4*j], l_h[2+4*j], l_h[3+4*j]
            
        l_h  = l_cont[i_line].split(); i_line += 1
        wyckoff_site_count[i] = l_h[0]

        for j in range(wyckoff_site_count[i]):
            l_h  = l_cont[i_line].split(); i_line += 1
            wyckoff_pos_count[j, i],  wyckoff_mult[j, i], \
                wyckoff_label[j, i] = l_h[0], l_h[1], l_h[2]

            for k in range(wyckoff_pos_count[j, i]):
                l_h  = l_cont[i_line].split(); i_line += 1
                wyckoff_bns_fract[:, k, j, i] = numpy.array(l_h[0:3],
                                                            dtype=int)
                wyckoff_bns_fract_denom[k, j, i] = l_h[3]
                wyckoff_bns_xyz[:, :, k, j, i] = numpy.transpose(
                    numpy.array(l_h[4:13], dtype=int).reshape(3,3))
                wyckoff_bns_mag[:, :, k, j, i] = numpy.transpose(
                    numpy.array(l_h[13:22], dtype=int).reshape(3,3))

        if magtype[i] == 4:
            l_h  = l_cont[i_line].split(); i_line += 1
            ops_count[i] = l_h[0]

            for j in range(ops_count[i]):
                if j%4 == 0:
                    l_h  = l_cont[i_line].split(); i_line += 1
                ops_og_point_op[j, i], ops_og_trans[0, j, i], \
                    ops_og_trans[1, j, i], ops_og_trans[2, j, i], \
                    ops_og_trans_denom[j, i], ops_og_timeinv[j, i] = \
                    l_h[0+6*j%4], l_h[1+6*j%4], l_h[2+6*j%4], l_h[3+6*j%4], \
                    l_h[4+6*j%4], l_h[5+6*j%4]
                
            l_h  = l_cont[i_line].split(); i_line += 1
            lattice_og_vectors_count[i] = l_h[0]

            l_h  = l_cont[i_line].split(); i_line += 1
            for j in range(lattice_og_vectors_count[i]):
                lattice_og_vectors[:, j, i] = numpy.array(l_h[0+4*j:3+4*j],
                                                          dtype=int)
                lattice_og_vectors_denom[j, i] = l_h[3+4*j]
                
            l_h  = l_cont[i_line].split(); i_line += 1
            wyckoff_site_count[i] = l_h[0]

            for j in range(wyckoff_site_count[i]):
                l_h  = l_cont[i_line].split(); i_line += 1
                wyckoff_pos_count[j, i], wyckoff_mult[j, i] = l_h[0], l_h[1]
                wyckoff_label[j, i] = l_h[2]
                
                for k in range(wyckoff_pos_count[j, i]):
                    l_h  = l_cont[i_line].split(); i_line += 1
                    wyckoff_og_fract[:, k, j, i] = numpy.array(l_h[0:3],
                                                               dtype=int)
                    wyckoff_og_fract_denom[k, j, i] = l_h[3]
                    wyckoff_og_xyz[:, :, k, j, i] = numpy.transpose(
                        numpy.array(l_h[4:4+9], dtype=int).reshape(3, 3))
                    wyckoff_og_mag[:, :, k, j, i] = numpy.transpose(
                        numpy.array(l_h[13:13+9], dtype=int).reshape(3, 3))

# read_magnetic_data()
# i = 5
# print("numerical label in BNS setting: ", nlabel_bns[i])
# print("jth part of nlabel_bns: ", nlabelparts_bns[:, i])
# print("group symbol: ", spacegroup_label_bns[i])
# print("numerical label in OG setting: ", nlabel_og[i])
# print("jth part of nlabel_og: ", nlabelparts_og[:, i])
# print("group symbol: ", spacegroup_label_og[i])
# print("type of magnetic space group (1-4): ", magtype[i])
# if magtype[i] == 4:
#     print("BNS-OG: 3x3 point operator part of transformation: ",
#           bnsog_point_op[:, :, i])
#     print("BNS-OG: translation part of transformation: ",
#           bnsog_origin[:, i])
#     print("BNS-OG: common denominator: ",
#           bnsog_origin_denom [i])
# print("number of point operators: ", ops_count[i])
# print("number of wyckoff sites: ", wyckoff_site_count[i])
# for j in range(wyckoff_site_count[i]):
#     print("number of positions in jth wyckoff site:", wyckoff_pos_count[j, i])
#     print("multiplicity for jth wyckoff site:", wyckoff_mult[j, i])
#     print("symbol (a,b,c,...,z,alpha) for jth wyckoff site:",
#           wyckoff_label[j, i])

#     print("kth component of fractional part of wyckoff position:",
#           wyckoff_bns_fract[:, :, j, i])
#     print("common denominator:", wyckoff_bns_fract_denom[:, j, i])


# print("number of lattice vectors defining the lattice: ",
#       lattice_bns_vectors_count[i])

# print("kth component of the jth lattice vector: ",
#       lattice_bns_vectors[:, :, i])
# print("common denominator: ", lattice_bns_vectors_denom[:, i])


# # for jth operator
# # ops_bns_point_op(j,i): point operator part
# ops_bns_point_op = numpy.zeros(shape = (96, magcount, ), dtype=int)

# # ops_bns_trans(k,j,i): kth component of translation part
# # ops_bns_trans_denom(j,i): common denominator
# ops_bns_trans = numpy.zeros(shape = (3, 96, magcount, ), dtype=int)
# ops_bns_trans_denom = numpy.zeros(shape = (96, magcount, ), dtype=int)

# # ops_bns_timeinv(j,i): 1=no time inversion, -1=time inversion
# ops_bns_timeinv = numpy.zeros(shape = (96, magcount, ), dtype=int)


# # wyckoff_bns_xyz(m,k,j,i): mth component to coeffcient of kth paramater (x,y,z)
# wyckoff_bns_xyz = numpy.zeros(shape = (3, 3, 96, 27, magcount, ), dtype=int)

# # wyckoff_bns_mag(m,k,j,i): mth component to coeffcient of kth magnetic
# # paramater (mx,my,mz)
# wyckoff_bns_mag = numpy.zeros(shape = (3, 3, 96, 27, magcount, ), dtype=int)


# # for OG setting (for type-4 groups)
# # lattice_og_vectors_count(i): number of lattice vectors defining the lattice
# lattice_og_vectors_count = numpy.zeros(shape = (magcount, ), dtype=int)

# # lattice_og_vectors(k,j,i): kth component of the jth lattice vector
# # lattice_og_vectors_denom(j,i): common denominator
# lattice_og_vectors = numpy.zeros(shape = (3, 6, magcount, ), dtype=int)
# lattice_og_vectors_denom = numpy.zeros(shape = (6, magcount, ), dtype=int)

# # for jth operator
# # ops_og_point_op(j,i): point operator part
# ops_og_point_op = numpy.zeros(shape = (96, magcount, ), dtype=int)

# # ops_og_trans(k,j,i): kth component of translation part
# # ops_og_trans_denom(j,i): common denominator
# ops_og_trans = numpy.zeros(shape = (3, 96, magcount, ), dtype=int)
# ops_og_trans_denom = numpy.zeros(shape = (96, magcount, ), dtype=int)

# # ops_og_timeinv(j,i): 1=no time inversion, -1=time inversion
# ops_og_timeinv = numpy.zeros(shape = (96, magcount, ), dtype=int)

# # for jth wyckoff site
# # wyckoff_og_fract(k,j,i): kth component of fractional part of wyckoff position
# # wyckoff_og_fract_denom(j,i): common denominator
# wyckoff_og_fract = numpy.zeros(shape = (3, 96, 27, magcount, ), dtype=int)
# wyckoff_og_fract_denom = numpy.zeros(shape = (96, 27, magcount, ), dtype=int)

# # wyckoff_og_xyz(m,k,j,i): mth component to coeffcient of kth paramater (x,y,z)
# wyckoff_og_xyz = numpy.zeros(shape = (3, 3, 96, 27, magcount, ), dtype=int)

# # wyckoff_og_mag(m,k,j,i): mth component to coeffcient of kth magnetic
# # paramater (mx,my,mz)
# wyckoff_og_mag = numpy.zeros(shape = (3, 3, 96, 27, magcount, ), dtype=int)

