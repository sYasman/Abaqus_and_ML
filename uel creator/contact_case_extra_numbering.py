### If additional contacts are created, after uel_prep this file must be interpreted
### Reading tset files create new element numbering
### TODO: must be implemented into uel_prep 

incr = 124308 # basically number of user elements in model

resp = input('Did you update the incr ??')

if resp == 'y':
    print(f'incr: {incr}')
    updated_file = 'new_001.inp'
    f_out = open(updated_file, 'w')
    with open('tset_001.inp', 'r') as f_in:
        r = 0
        start_found = False
        for r, lines in enumerate(f_in):
            r += 1
            if start_found:
                if '*' in lines:
                    f_out.write(f'\n{lines}\n')
                    start_found = False       
                else:
                    aa = lines.strip().split(',')
                    i = 0
                    for a in aa:
                        if a != '':
                            i += 1
                            # print(int(a)+incr,end=',')
                            f_out.write(f'{int(a)+incr},')
                            if i > 5:
                                f_out.write('\n')
                                i = 0
                    # break
            else:
                # print(lines, end='')
                f_out.write(f'\n{lines}\n')
                
            if '*Elset' in lines and 'MATRIX_WITH_HOLES-1' in lines:
                start_row = r
                start_found =True
                # print(lines, r)
    f_out.close()

    print(f'Fine Updated: {updated_file}')
else:
    print('Update incr...')
	
import sys
output=""
with open(updated_file) as f:
    for line in f:
        if not line.isspace():
            output+=line
            
f= open(updated_file,"w")
f.write(output)	