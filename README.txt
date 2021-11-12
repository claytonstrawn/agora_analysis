# agora_analysis

To install this file, clone from the github account and then use pip.

git clone https://github.com/claytonstrawn/agora_analysis.git
cd agora_analysis
pip install -e .

Create snapshots for different codes by using 

snap = agora_analysis.AgoraSnapshot(name, redshift)

where "name" is of the form "AGORA_<code>_<simulation>", where <code> is one of 'art','enzo','ramses','gadget','gear','gizmo','changa', and <simulation> is one of 'CR' (CosmoRun), 'C1' (Cal1), 'C2' (Cal2), or 'C3' (Cal3). [currently only 'CR' and 'C1' are implemented].

Redshift is any float from z=11.0 to z=2.0, and it will look up the closest snapshot to that value.

Once the snap is created, it can look up metadata like the file location of the snap, the location of the center of the most massive galaxy in kpc, and the virial radius, either specific to that code or the average virial radius between all codes in that simulation at that time. Ex:

snap = agora_analysis.AgoraSnapshot('AGORA_art_CR',4.0)
snap.center
>>>unyt_array([8435.8744002 , 8880.40809419, 8600.03624835], 'kpc')
snap.Rvir
>>>unyt_quantity(21.1835916, 'kpc')
snap.lookup_snap_path()
>>>'/project/projectdirs/agora/paper_CGM/ThinTimeSteps/ART-I/Cosmo_v14_final/10MpcBox_csf512_01154.d'

To analyze further using YT and compare multiple codes at once, you can load the snap:

snap.load_snapshot()

Which will automatically look up the file location and load the snapshot. All analysis scripts for Paper IV, and hopefully we will be able to backport all scripts from Paper III, II, and I, will be available in image_scripts for generating the figures from those papers directly, or analysis_scripts for generating the numbers which are used in more quantitative images. 

Field_setup contains methods for creating uniform yt fields between codes: using these each code can have a comparable field ('agora_stars','agora_mass'), ('gas','agora_density'), etc. instead of each code having a separately named set of fields. [this part not fully implemented/explained yet, will edit with development.]