data = open("InputParameters.txt", "r")
lines = data.readlines()
data.close()

lst = []
for line in lines:
    a = line.strip("\n")    #removes the \n
    b = a.replace("</br>", " ")     #replaces all the /br with blank spaces so that I can take pretty much every word separetely
    c = b.split()       #separates all the words
    for word in c:
        lst.append(word)

prop = lst[::4]

symbol = lst[1::4]

value = lst[2::4]

unit = lst[3::4]
