// Search for the appropriate blockMesh dictionary
const word dictName("blockMeshDict");

autoPtr<IOdictionary> meshDictPtr;

{
    fileName dictPath;

    if (args.optionFound("dict"))
    {
        // Dictionary specified on the command-line ...

        if (isDir(dictPath))
        {
            dictPath = dictPath/dictName;
        }
    }
    else if
    (
        exists
        (
            runTime.path()/runTime.constant()
            /regionPath/polyMesh::meshSubDir/dictName
        )
    )
    {
        // Dictionary present in constant polyMesh directory (old-style)

        dictPath = runTime.constant()/regionPath/polyMesh::meshSubDir/dictName;

        // Warn that constant/polyMesh/blockMeshDict was used
        // instead of system/blockMeshDict
        WarningInFunction
            << "Using the old blockMeshDict location: "
            << dictPath << nl
            << "    instead of the default location:  "
            << runTime.system()/regionPath/dictName << nl
            << endl;
    }
    else
    {
        // Assume dictionary is to be found in the system directory

        dictPath = runTime.system()/regionPath/dictName;
    }

    IOobject meshDictIO
    (
        dictPath,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (!meshDictIO.headerOk())
    {
        FatalErrorInFunction
            << "Cannot open blockMeshDict "
            << meshDictIO.objectPath() << nl
            << exit(FatalError);
    }

    Info<< nl << "Creating block mesh from\n    "
        << meshDictIO.objectPath() << nl << endl;

    meshDictPtr = autoPtr<IOdictionary>(new IOdictionary(meshDictIO));
}

const IOdictionary& meshDict = meshDictPtr;
