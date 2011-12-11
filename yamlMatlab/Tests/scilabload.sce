jimport java.io.InputStream;
jimport java.io.File;
res=javaclasspath();
javaclasspath('/home/jirka/Software/yamlmatlab/trunk/external-packages/snakeyaml/snakeyaml-1.8.jar');

jimport('org.yaml.snakeyaml.Yaml')

yamlreader = Yaml.new();
u = mopen('/home/jirka/SVN/EnergoCentrum/MPC_RT/Dejvice/Definition/OptimalControlProblem.yml','r');
yml = mgetl(u);
mclose(u);
jymlobj = yamlreader.load(yml);

f = File.new("home/jirka/SVN/EnergoCentrum/MPC_RT/Dejvice/Definition/OptimalControlProblem.yml");
 input = FileInputStream.new(f);

