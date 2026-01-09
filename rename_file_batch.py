import os

def rename_file_batch(folder_path,change_from,change_to):

    # Set the directory containing your .lnk files
    # folder_path = r"C:\Users\hanshil\OneDrive - NOC\Documents\Ocean_Optics\SINK_Data\SINK2025\LISST\20250601_reprocess"

    for filename in os.listdir(folder_path):
        if change_from in filename:
            original_path = os.path.join(folder_path, filename)
            new_filename = filename.replace(change_from, change_to)
            new_path = os.path.join(folder_path, new_filename)

            os.rename(original_path, new_path)
            print(f"Renamed: {filename} -> {new_filename}")
