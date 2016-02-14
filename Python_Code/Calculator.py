
def add(num1, num2):
	return num1 + num2
	
def subtract(num1, num2):
	return num1 - num2
	
def multiply(num1, num2):
	return num1 * num2
	
def divide(num1, num2):
	return num1 / num2
	
def main():
	opperation = input("What do you want to do? (+,-,*,/): ")
	if ((opperation != '+')&(opperation != '-')&(opperation != '*')&(opperation != '/')):
		print("Please enter one of the values in brackets ")
		main()
	else:
		num1 = float(input("Enter the first number: "))
		num2 = float(input("Enter the second number: "))
		if(opperation == '+'):
			print(add(num1, num2))
		elif(opperation == '-'):
			print(subtract(num1, num2))
		elif(opperation == '*'):
			print(multiply(num1, num2))
		elif(opperation == '/'):
			print(divide(num1, num2))
	
main()