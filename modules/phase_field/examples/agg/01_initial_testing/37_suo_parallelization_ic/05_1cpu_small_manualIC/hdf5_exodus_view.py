import h5py

def print_name(name, obj):
    print(name)

# Open the HDF5 file in read mode
with h5py.File('05_1cpu_small_manualIC_hdf5_exodus.e', 'r') as f:
    # .visititems visits every group and dataset in the file
    f.visititems(print_name)
