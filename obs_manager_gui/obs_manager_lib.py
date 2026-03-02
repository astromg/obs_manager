class MiscLine:
    def __init__(self, raw_line, line_index):
        self.ob_line = False
        self.raw_line = raw_line
        self.line_index = line_index




class ObjectLine(MiscLine):
    def __init__(self, raw_line, line_index, columns):
        self.ob_line = True
        self.raw_line = raw_line
        self.line_index = line_index

        # ta czesc jest tymczasowa, do uzycia parsera

        tmp = {}
        line = raw_line

        tmp["name"] = raw_line.split()[0]
        tmp["ra"] = raw_line.split()[1]
        tmp["dec"] = raw_line.split()[2]

        tmp["other"] = []
        tmp["editted"] = []

        for x in columns:
                if x + "=" in line and "comment" not in x:
                    tmp[x] = line.split(x + "=")[1].split()[0]
                    line = line.replace(f'{x}={tmp[x]}', "")

                if "comment=" in line:
                    comment = line.split("comment=")[1]
                    if len(comment) > 0:
                        if "\"" == comment[0]:
                            comment = comment.split("\"")[1]
                            line = line.replace(f'comment=\"{comment}\"', "")
                        elif "(" == comment[0]:
                            comment = comment.split("(")[1].split(")")[0]
                            line = line.replace(f'comment=({comment})', "")
                        else:
                            comment = comment.split()[0]
                            line = line.replace(f'comment={comment}', "")
                    tmp["comment"] = comment
                tmp["other"] = line.strip()

        self.ob = tmp



class MasterList:
    def __init__(self):
        self.lines=[]
    def add_line(self, line, columns, i = None):
        ob_line = False
        if len(line) > 0:
            if len(line.split()) > 0:
                if "#" not in line.split()[0]:
                    ob_line = True

        if i is not None:
            if ob_line:
                tmp = ObjectLine(line, i, columns)
            else:
                tmp = MiscLine(line, i)
            self.lines.insert(i, tmp)
            for a in self.lines:
                if a.line_index > i:
                    a.line_index += 1

        else:
            i = len(self.lines)
            if ob_line:
                tmp = ObjectLine(line, i, columns)
            else:
                tmp = MiscLine(line, i)
            self.lines.append(tmp)