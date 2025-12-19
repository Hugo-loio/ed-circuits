using ProfileView
using FileIO

profile_name = "haar_circuit.jlprof"

data = load(profile_name)
ProfileView.view(data[1], lidict=data[2])
