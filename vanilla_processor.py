# simple_processor.py

from coffea.processor import ProcessorABC
import awkward as ak

class PrintFieldsProcessor(ProcessorABC):
    def process(self, events):
        print("Fields in this event file:")
        print(events.fields)
        return {}

    def postprocess(self, accumulator):
        return accumulator
