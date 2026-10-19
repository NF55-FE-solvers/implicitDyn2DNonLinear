# this code generates the list of nodes given in Abaqus .inp file implicitly
# by implicitly I mean: start, end, delta. This is typically specified with 
# keyword "generate"

# Insert your data
start = int(input("Insert the STARTing nodeID: "))
end   = int(input("Insert the ENDing nodeID: "))
delta = int(input("Insert the DELTA: "))

# Printing the list
val = start
print(f"{val}", end=",")
while(val < end):
    val = val + delta
    print(f"{val}", end=",")

