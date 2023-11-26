
filename = 'SmartField4gau.gjf'
ff_string=[]
with open(filename, 'r') as file:
        file_content = file.read()
        split_content = file_content.split('!Master function')
        if len(split_content) > 1:
            # Storing the text after encountering the '!Master function'
            text_after_master_function = split_content[1].strip()
            ff_string.append(text_after_master_function)
print(ff_string)
