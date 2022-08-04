import time
import socket

with open(snakemake.output[0], "w") as f:
    msg = f"Hello from {socket.gethostname()}!\n"
    f.write(msg)
    time.sleep(2)
