/*---------------------------------------------------------------------------*\
|                       _    _  _     ___                       | The         |
|     _____      ____ _| | _| || |   / __\__   __ _ _ __ ___    | Swiss       |
|    / __\ \ /\ / / _` | |/ / || |_ / _\/ _ \ / _` | '_ ` _ \   | Army        |
|    \__ \\ V  V / (_| |   <|__   _/ / | (_) | (_| | | | | | |  | Knife       |
|    |___/ \_/\_/ \__,_|_|\_\  |_| \/   \___/ \__,_|_| |_| |_|  | For         |
|                                                               | OpenFOAM    |
-------------------------------------------------------------------------------
License
    This file is part of swak4Foam.

    swak4Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    swak4Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with swak4Foam; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Contributors/Copyright:
    2008-2011, 2013, 2015-2016, 2018 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>

 SWAK Revision: $Id$
\*---------------------------------------------------------------------------*/

#include "writeAdditionalFieldsFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

#include "polyMesh.H"
#include "IOmanip.H"
#include "swakTime.H"
// #include "IOobjectOption.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(writeAdditionalFieldsFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeAdditionalFieldsFunctionObject,
        dictionary
    );

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

writeAdditionalFieldsFunctionObject::writeAdditionalFieldsFunctionObject
(
    const word &name,
    const Time& t,
    const dictionary& dict
)
:
    writeFieldsGeneralFunctionObject(
        name,t,dict,false
    ),
    no_write_(
        dict.lookupOrDefault<bool>("no_write", false)
    )
{
}

bool writeAdditionalFieldsFunctionObject::start()
{
    writeFieldsGeneralFunctionObject::start();

    label found = 0;
    label changed = 0;

    // IOobjectOption::writeOption expected_setting = no_write_ ? IOobject::NO_WRITE : IOobject::AUTO_WRITE;
    IOobject::writeOption expected_setting = no_write_ ? IOobject::NO_WRITE : IOobject::AUTO_WRITE;

    forAll(fieldNames(), i) {
        const word &name=fieldNames()[i];
        if(obr_.foundObject<IOobject>(name)) {
            found ++;
            IOobject &obj=const_cast<IOobject&>(
                obr_.lookupObject<IOobject>(name)
            );
            if(obj.writeOpt() != expected_setting) {
                changed ++;
                obj.writeOpt() = expected_setting;
            }
        }
    }

    Info << "Additional fields " << fieldNames() << " will "
        << (no_write_ ? "not " : "") << "be written. "
        << found << " of these are already present and " << changed
        << " were changed to "
        << (no_write_ ? "NO_WRITE" : "AUTO_WRITE")
        << endl;

    if(outputControlMode()==ocmStartup) {
        writeSimple();
    }

    return true;
}

// bool writeAdditionalFieldsFunctionObject::outputTime()
// {
//     return (
//         time().outputTime()
//         &&
//         time().time().value()>=after());
// }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

} // namespace Foam

// ************************************************************************* //
