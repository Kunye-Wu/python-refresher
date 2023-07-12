class BankAccount:
    def __init__(self, balance, name, account_number):
        self.balance = balance
        self.name = name
        self.account_number = account_number

    def withdraw(self, amount):
        if amount > self.balance:
            print("Insufficient funds. Cannot withdraw.")
        else:
            self.balance -= amount
            print(f"Withdrawal successful. New balance: {self.balance}")

    def deposit(self, amount):
        self.balance += amount
        print(f"Deposit successful. New balance: {self.balance}")

    def print_balance(self):
        print(
            f"Account: {self.account_number}\nName: {self.name}\nBalance: {self.balance}"
        )


if __name__ == "__main__":
    # Usage example:
    p1 = BankAccount(1000, "John Chen", "123456")
    p1.print_balance()  # Print initial balance
    p1.withdraw(800)  # Withdraw 800
    p1.deposit(1500)  # Deposit 1500
    p1.print_balance()  # Print updated balance
