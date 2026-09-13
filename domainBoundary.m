common_data_dir  = getenv("COMMON_DATA");
mask_folder = common_data_dir + "drainage-basins/MEaSURE-Antarctic-boundaries/";

imbie_path = mask_folder + "Basins_IMBIE_Antarctica_v02.dbf";
info = shaperead(convertStringsToChars(imbie_path));
