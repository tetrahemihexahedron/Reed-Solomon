import RS
import FiniteFields as FF
from Tkinter import *

class DecoderWindow(Tk):
    def __init__(self, parent):
        Tk.__init__(self, parent)
        self.parent = parent
        self.title("Reed-Solomon Decoding")
        
        self.initialize()

    def initialize(self):
        instructionsFrame = Frame(self)
        instructions = Label(instructionsFrame,
                             text = "Input the received transmission.",
                             anchor = NW)
        instructions.pack(fill = X, padx = 5, pady = 5)
        instructionsFrame.pack(fill = BOTH)

        transmissionFrame = Frame(self)
        firstEntryLabel = Label(transmissionFrame, text = "m(0) = ")
        firstEntryLabel.grid(row = 0, column = 0, sticky = E)
        self.firstEntry = Entry(transmissionFrame, width = 5)
        self.firstEntry.grid(row = 0, column = 1)

        secondEntryLabel = Label(transmissionFrame, text = "m(1) = ")
        secondEntryLabel.grid(row = 0, column = 2, sticky = E)
        self.secondEntry = Entry(transmissionFrame, width = 5)
        self.secondEntry.grid(row = 0, column = 3)

        thirdEntryLabel = Label(transmissionFrame, text = "m(t) = ")
        thirdEntryLabel.grid(row = 0, column = 4, sticky = E)
        self.thirdEntry = Entry(transmissionFrame, width = 5)
        self.thirdEntry.grid(row = 0, column = 5)

        fourthEntryLabel = Label(transmissionFrame, text = "m(t^2) = ")
        fourthEntryLabel.grid(row = 0, column = 6, sticky = E)
        self.fourthEntry = Entry(transmissionFrame, width = 5)
        self.fourthEntry.grid(row = 0, column = 7)

        fifthEntryLabel = Label(transmissionFrame, text = "m(t^3) = ")
        fifthEntryLabel.grid(row = 1, column = 0, sticky = E)
        self.fifthEntry = Entry(transmissionFrame, width = 5)
        self.fifthEntry.grid(row = 1, column = 1)

        sixthEntryLabel = Label(transmissionFrame, text = "m(t^4) = ")
        sixthEntryLabel.grid(row = 1, column = 2, sticky = E)
        self.sixthEntry = Entry(transmissionFrame, width = 5)
        self.sixthEntry.grid(row = 1, column = 3)

        seventhEntryLabel = Label(transmissionFrame, text = "m(t^5) = ")
        seventhEntryLabel.grid(row = 1, column = 4, sticky = E)
        self.seventhEntry = Entry(transmissionFrame, width = 5)
        self.seventhEntry.grid(row = 1, column = 5)

        eighthEntryLabel = Label(transmissionFrame, text = "m(t^6) = ")
        eighthEntryLabel.grid(row = 1, column = 6, sticky = E)
        self.eighthEntry = Entry(transmissionFrame, width = 5)
        self.eighthEntry.grid(row = 1, column = 7)
        
        transmissionFrame.pack(ipadx = 20)

        decodeButtonFrame = Frame(self)
        decodeButton = Button(decodeButtonFrame, text = "Decode", command = self.decode)
        decodeButton.pack(pady = 5, padx = 5)
        decodeButtonFrame.pack(fill = BOTH)

        resultsFrame = Frame(self)
        self.resultsText = StringVar()
        results = Label(resultsFrame, textvariable = self.resultsText, anchor = NW)
        results.pack(fill = X, pady = 5, padx = 5)
        resultsFrame.pack(fill = BOTH)

    def decode(self):
        T = self.getTransmission()
        results = RS.RSDecode(tuple(T), 4)
        self.resultsText.set(str(results))

    def getTransmission(self):
        e1 = self.firstEntry.get()
        e2 = self.secondEntry.get()
        e3 = self.thirdEntry.get()
        e4 = self.fourthEntry.get()
        e5 = self.fifthEntry.get()
        e6 = self.sixthEntry.get()
        e7 = self.seventhEntry.get()
        e8 = self.eighthEntry.get()
        transmission = []
        for e in (e1, e2, e3, e4, e5, e6, e7, e8):
            if e == "1":
                transmission += [FF.FFieldElt(FF.FiniteField((1, 1, 0, 1)),
                                              [1, 0, 0])]
            elif e == "0":
                transmission += [FF.FFieldElt(FF.FiniteField((1, 1, 0, 1)),
                                             [0, 0, 0])]
            else:
                poly = []
                for c in e:
                    poly += [int(c)]
                poly.reverse()
                transmission += [FF.FFieldElt(FF.FiniteField((1, 1, 0, 1)),
                                         poly)]
        return transmission
        
    
if __name__ == "__main__":
    decoder = DecoderWindow(None)
    decoder.mainloop()
