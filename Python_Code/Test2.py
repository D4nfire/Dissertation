
#This is a comment.

def helloworld(myString):
	print(myString)
	myName = input("What is your name? ")	
	print("Hello " + myName)

	#You can use or, or |
	
	if((myName == 'dan')|(myName == 'Dan')):
		print("You are awesome! ")
	elif(myName == 'lee' or myName == 'Lee'):
		print("You are da bomb! ")
	else:
		print("You are alright ")
	
helloWorld("Hello world ")