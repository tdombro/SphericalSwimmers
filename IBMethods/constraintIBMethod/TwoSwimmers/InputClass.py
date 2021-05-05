import sys, os
from utilities import remove_comment



class InputDatabase:
    "store input parameter from file"
    def __init__(self, filename):
        self.d_filename = filename
        self.d_list = {}
        self.d_inner_database = {}
        if filename != '':
            self.setDataFile(filename)

    def getDatabase(self, subkey):
        return self.d_inner_database[subkey]

    def setDataString(self, raw_lines):
        #print(raw_lines)
        i = 0
        #print(' len(raw_lines) = ', len(raw_lines))
        while i < len(raw_lines):
            proc_line = remove_comment(raw_lines[i])
            # splitting to words
            words = proc_line.split()
            lw = len(words)
            i += 1
            if lw == 0: continue
            if lw >= 3:
                # set key = value line
                self.d_list[words[0]] = words[2:]
                # sub-key line
            elif lw == 1:
                subkey = words[0]
                # if the next line not starting with {, then invalid file
                # note i already increases by 1
                if raw_lines[i][0] != '{':
                    # print('level = ', self.d_level)
                    print(i, raw_lines[i], " invalid input data")
                    quit()
                # copying from i + 1 until } met
                sub_lines = []
                i += 1
                while raw_lines[i][0] != '}':
                    sub_lines.append(raw_lines[i])
                    i += 1
                i += 1
                # make a child database
                # print(sub_lines)
                if len(sub_lines) != 0:
                    inner_database = InputDatabase('')
                    inner_database.setDataString(sub_lines)
                    self.d_inner_database[subkey] = inner_database

    def setDataFile(self, filename):
        # read input from file in the current dir
        current_dir = os.getcwd() + '/'
        print(current_dir + self.d_filename)
        input_dat = open(current_dir + self.d_filename, 'r')
        raw_lines = input_dat.readlines()
        self.setDataString(raw_lines)

        # process each line at a time

    def getBool(self, key):
        v = self.getString(key)
        value = v.lower()
        if value == 'y' or value == 'yes':
            return True
        else:
            return False

    def getDouble(self, key):
        value = self.d_list[key]
        return float(value[0])

    def getInt(self, key):
        value = self.d_list[key]
        return int(value[0])

    def getDoubleArray(self, key):
        value = self.d_list[key]
        ret = []
        for x in value: ret.append(float(x))
        return ret

    def getString(self, key):
        s = ' '.join(self.d_list[key])
        if s[0] == '"' or s[0] == "'":
            return s[1:-1]
        else:
            return s



    def print(self):
        for key in self.d_list:
            print(key, " = ", self.d_list[key])


#p = InputDatabase('input.spherobot')
#a = p.getDatabase('genSwimmers').getString('run')
#print(a)
#print(p.getString('test'))
#print(p.getString('post_dir_name'))


