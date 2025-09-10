# Define bit meanings
bit_meanings = {
    0: "isPrompt",
    1: "isDecayedLeptonHadron",
    2: "isTauDecayProduct",
    3: "isPromptTauDecayProduct",
    4: "isDirectTauDecayProduct",
    5: "isDirectPromptTauDecayProduct",
    6: "isDirectHadronDecayProduct",
    7: "isHardProcess",
    8: "fromHardProcess",
    9: "isHardProcessTauDecayProduct",
    10: "isDirectHardProcessTauDecayProduct",
    11: "fromHardProcessBeforeFSR",
    12: "isFirstCopy",
    13: "isLastCopy",
    14: "isLastCopyBeforeFSR",
}

# Example input list of decimal numbers
decimal_numbers =[ 8449, 12288, 12289, 12308, 12348, 12352, 12356, 12364, 12673, 22913]

# Write results to file
with open("gen_status_flags.txt", "w") as f:
    for num in decimal_numbers:
        binary_str = format(num, "015b")  # 15-bit binary
        f.write(f"Decimal: {num}\n")
        f.write(f"Binary : {binary_str}\n")
        f.write("Properties:\n")
        
        for bit, meaning in bit_meanings.items():
            if num & (1 << bit):
                f.write(f"  - {meaning}\n")
        f.write("\n")

