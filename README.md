On Ubuntu 24.04 LTS:

```bash
$ sudo apt install make g++ libglfw3-dev libglew-dev libglm-dev libluajit-5.1-dev libgmp-dev libsqlite3-dev
$ make
$ make run
```

For clang:
```bash
$ make barbarize
$ make all BARBARIZED=true CXX=clang++
$ make run
```

Don’t forget to copy the configuration file:
```bash
$ cp config.lua.example config.lua
```

## License

For any copyright year range specified as YYYY–ZZZZ in this package note
that the range specifies every single year in that closed interval.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
