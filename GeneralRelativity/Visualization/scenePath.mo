within GeneralRelativity.Visualization;
impure function scenePath "Full path of a scene file: chooses the output folder and creates it if needed"
  input String fileName "Scene file name, for example schwarzschild_scene.ndjson";
  input String directory = "" "Output folder; empty selects the default described below";
  output String path "Folder plus file name, or the bare file name for the working directory";
protected
  String chosen;
  String overrideDirectory;
  Boolean overrideSet;
  String home;
  Boolean homeSet;
algorithm
  chosen := directory;
  if chosen == "" then
    (overrideDirectory, overrideSet) := Modelica.Utilities.System.getEnvironmentVariable("GR_SCENE_DIR", true);
    if overrideSet and overrideDirectory <> "" then
      chosen := overrideDirectory;
    else
      (home, homeSet) := Modelica.Utilities.System.getEnvironmentVariable("HOME", true);
      if not homeSet or home == "" then
        (home, homeSet) := Modelica.Utilities.System.getEnvironmentVariable("USERPROFILE", true);
      end if;
      // Only use Documents when it already exists, so headless runners with an
      // unwritable or bare home directory fall back to the working directory.
      if homeSet and home <> "" and Modelica.Utilities.Files.exist(home + "/Documents") then
        chosen := home + "/Documents/Visualization";
      end if;
    end if;
  end if;
  if chosen == "" then
    path := fileName;
  else
    Modelica.Utilities.Files.createDirectory(chosen);
    path := chosen + "/" + fileName;
  end if;
  annotation(Documentation(info="<html>
<p>The folder is chosen in this order:</p>
<ol>
<li><code>directory</code>, when it is not empty;</li>
<li>the environment variable <code>GR_SCENE_DIR</code>, when it is set;</li>
<li><code>~/Documents/Visualization</code> (from <code>HOME</code>, or <code>USERPROFILE</code> on
Windows), when <code>~/Documents</code> exists;</li>
<li>otherwise the simulation working directory.</li>
</ol>
<p>Missing folders are created. Use forward slashes; they work on every platform.</p>
</html>"));
end scenePath;
